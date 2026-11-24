"""Sequential and process-parallel pruning for OPS route search."""

from __future__ import annotations

from array import array
from concurrent.futures import ProcessPoolExecutor, wait
from contextlib import contextmanager
import math
import mmap
from multiprocessing import get_context
import os
import shutil
import tempfile
from types import SimpleNamespace
from typing import TYPE_CHECKING

from OPS.src.write_progress import DiskRouteSet, RouteChunkWriter, _iter_disk_routes

if TYPE_CHECKING:
    from OPS.src.basic_class import remezData


ROUTE_CHUNK_SIZE = 100_000
PARALLEL_MIN_ROUTES = 10_000

# BEST vs candidates
def prune_with_BEST(
    routes: list[remezData] | DiskRouteSet,
    BEST: remezData,
    options: SimpleNamespace,
    *,
    depth_margin: int = 2,
) -> tuple[list[remezData] | DiskRouteSet, int]:
    orig_count = len(routes)

    if options.criteria == "time":
        keep_route = lambda r: r.total_time <= BEST.total_time
    elif options.criteria == "depth":
        keep_route = lambda r: r.total_CalData.depth <= BEST.total_CalData.depth - depth_margin
    else:
        raise ValueError(f"Unsupported criteria: {options.criteria}")

    if isinstance(routes, DiskRouteSet):
        # DiskRouteSet is intentionally not directly iterable. Filter each
        # persisted chunk so large searches stay within the memory bound.
        parent_dir = os.path.dirname(routes.dir_path)
        pruned_dir = tempfile.mkdtemp(prefix="pruned_", dir=parent_dir)
        writer = RouteChunkWriter(pruned_dir, ROUTE_CHUNK_SIZE)

        try:
            for _, chunk_routes in _iter_disk_routes(routes):
                for route in chunk_routes:
                    if keep_route(route):
                        writer.add(route)
            filtered_routes = writer.close()
        except Exception:
            shutil.rmtree(pruned_dir, ignore_errors=True)
            raise

        routes.cleanup()
        routes = filtered_routes
    else:
        routes = [route for route in routes if keep_route(route)]

    rem_count = orig_count - len(routes)

    return routes, rem_count

# candidates vs candidates
def prune_with_candidates(
    routes: list[remezData],
    criteria: str,
    prune_range: tuple[str, ...],
    *,
    executor: ProcessPoolExecutor | None = None,
    workers: int = 1,
    parallel_min_routes: int = PARALLEL_MIN_ROUTES,
    pruning_type: str = "prune",
) -> tuple[list[remezData], int]:
    if criteria != "depth":
        raise ValueError(f"Unsupported criteria: {criteria}")
    if pruning_type not in ("prune", "lex"):
        raise ValueError("pruning_type must be one of: prune, lex")
    if workers < 1 or parallel_min_routes < 0:
        raise ValueError("workers must be positive and parallel_min_routes nonnegative")
    if workers == 1 or len(routes) < max(2, parallel_min_routes):
        return _filter_depth_dominated(routes, prune_range, pruning_type)
    return prune_with_candidates_parallel(
        routes, criteria, prune_range, executor=executor, workers=workers,
        pruning_type=pruning_type,
    )

class _DominanceNode:
    __slots__ = ("minimum", "maximum", "left", "right", "points", "size")

    LEAF_SIZE = 16

    def __init__(
        self,
        points: list[tuple[int, tuple[float, ...]]],
        scales: tuple[float, ...] | None = None,
    ):
        dimension_count = len(points[0][1])
        self.minimum = tuple(min(point[dimension] for _, point in points) for dimension in range(dimension_count))
        self.maximum = tuple(max(point[dimension] for _, point in points) for dimension in range(dimension_count))
        self.left = None
        self.right = None
        self.points = None
        self.size = len(points)

        if len(points) <= self.LEAF_SIZE:
            self.points = points
            return

        if scales is None:
            scales = tuple(
                maximum - minimum if maximum > minimum else 1.0
                for minimum, maximum in zip(self.minimum, self.maximum)
            )

        ranges = [
            (self.maximum[dimension] - self.minimum[dimension]) / scales[dimension]
            for dimension in range(dimension_count)
        ]
        split_dimension = max(range(dimension_count), key=ranges.__getitem__)
        points.sort(key=lambda item: item[1][split_dimension])
        middle = len(points) // 2
        self.left = _DominanceNode(points[:middle], scales)
        self.right = _DominanceNode(points[middle:], scales)

    def _all_points_dominate(
        self,
        target: tuple[float, ...],
        exclude_index: int | None,
        require_point_strict: bool,
        gap_dimension: int | None,
        minimum_gap: float,
    ) -> bool:
        if exclude_index is not None and self.size == 1:
            return False

        for dimension, (maximum, limit) in enumerate(zip(self.maximum, target)):
            if dimension == gap_dimension:
                if limit - maximum <= minimum_gap:
                    return False
            elif maximum > limit:
                return False

        return not require_point_strict or self.minimum != target

    def has_dominator(
        self,
        target: tuple[float, ...],
        exclude_index: int | None,
        require_point_strict: bool,
        gap_dimension: int | None,
        minimum_gap: float,
    ) -> bool:
        for dimension, (minimum, limit) in enumerate(zip(self.minimum, target)):
            if dimension == gap_dimension:
                if limit - minimum <= minimum_gap:
                    return False
            elif minimum > limit:
                return False

        if self._all_points_dominate(
            target,
            exclude_index,
            require_point_strict,
            gap_dimension,
            minimum_gap,
        ):
            return True

        if self.points is not None:
            for index, point in self.points:
                if index == exclude_index:
                    continue
                if all(value <= limit for value, limit in zip(point, target)):
                    if gap_dimension is not None and target[gap_dimension] - point[gap_dimension] <= minimum_gap:
                        continue
                    if not require_point_strict or point != target:
                        return True
            return False

        children = (self.left, self.right)
        if self.right.minimum < self.left.minimum:
            children = (self.right, self.left)

        return any(
            child.has_dominator(
                target,
                exclude_index,
                require_point_strict,
                gap_dimension,
                minimum_gap,
            )
            for child in children
        )


class _DominanceBucketIndex:
    def __init__(
        self,
        entries: list[tuple[int, tuple[float, ...]]],
    ):
        self.root = _DominanceNode(entries)

    def has_dominator(
        self,
        target: tuple[float, ...],
        *,
        exclude_index: int | None = None,
        require_point_strict: bool = False,
        gap_dimension: int | None = None,
        minimum_gap: float = 0.0,
    ) -> bool:
        return self.root.has_dominator(
            target,
            exclude_index,
            require_point_strict,
            gap_dimension,
            minimum_gap,
        )

# depth criteria: DCPA buckets with a 4-D precision-containment index
def _prepare_depth_index(
    routes: list[remezData],
    prune_range: tuple[str, ...],
    pruning_type: str = "prune",
):
    def get_value(route: remezData):
        (lb0, ub0), (lb1, ub1) = route.interval_log[-1]
        cal_data = route.total_CalData
        optional_costs = {
            "c": cal_data.cmult,
            "p": cal_data.pmult,
            "a": cal_data.add_count,
        }
        cost_key = (
            cal_data.depth,
            *(optional_costs[attribute] for attribute in prune_range),
        )

        # All four coordinates are minimized. Component-wise comparison is
        # therefore equivalent to interval containment in both output ranges.
        precision = (-lb0, ub0, -lb1, ub1)
        return precision, cost_key

    values = [get_value(route) for route in routes]

    bucket_entries: dict[
        tuple[int, ...],
        list[tuple[int, tuple[float, float, float, float]]],
    ] = {}
    for index, (precision, cost_key) in enumerate(values):
        bucket_entries.setdefault(cost_key, []).append((index, precision))

    bucket_indices = {
        cost_key: _DominanceBucketIndex(entries)
        for cost_key, entries in bucket_entries.items()
    }

    cost_keys = list(bucket_entries)
    dominating_cost_keys: dict[tuple[int, ...], list[tuple[int, ...]]] = {}
    for target_cost in cost_keys:
        if pruning_type == "lex":
            # Keep the selected D, C, P, A order; interval containment is unchanged.
            compatible_costs = [source_cost for source_cost in cost_keys if source_cost <= target_cost]
        else:
            compatible_costs = [
                source_cost
                for source_cost in cost_keys
                if all(source <= target for source, target in zip(source_cost, target_cost))
            ]
        compatible_costs.sort(key=lambda cost: (sum(cost), cost))
        dominating_cost_keys[target_cost] = compatible_costs

    return values, bucket_indices, dominating_cost_keys


def _filter_depth_dominated(
    routes: list[remezData],
    prune_range: tuple[str, ...],
    pruning_type: str = "prune",
) -> tuple[list[remezData], int]:
    orig_count = len(routes)
    if orig_count < 2:
        return list(routes), 0

    values, bucket_indices, dominating_cost_keys = _prepare_depth_index(routes, prune_range, pruning_type)

    dominated = [False] * orig_count
    for target_index, (target_precision, target_cost) in enumerate(values):
        for source_cost in dominating_cost_keys[target_cost]:
            same_cost = source_cost == target_cost
            if bucket_indices[source_cost].has_dominator(
                target_precision,
                exclude_index=target_index if same_cost else None,
                require_point_strict=same_cost,
            ):
                dominated[target_index] = True
                break

    pruned_routes = [
        route for index, route in enumerate(routes)
        if not dominated[index]
    ]
    rem_count = orig_count - len(pruned_routes)

    return pruned_routes, rem_count


@contextmanager
def _shared_depth_index(values, bucket_indices, dominating_cost_keys):
    """Pack a read-only index in one anonymous Linux file shared by all workers."""
    if not hasattr(os, "memfd_create"):
        raise RuntimeError("Parallel OPS pruning requires Linux memfd_create and /proc")
    bounds = array("d")
    nodes = array("q")
    leaves = array("q")

    def flatten(node):
        index = len(nodes) // 5
        bounds.extend((*node.minimum, *node.maximum))
        nodes.extend((-1, -1, -1, -1, node.size))
        if node.points is not None:
            nodes[index * 5 + 2] = len(leaves)
            leaves.extend(point_index for point_index, _ in node.points)
            nodes[index * 5 + 3] = len(leaves)
        else:
            nodes[index * 5] = flatten(node.left)
            nodes[index * 5 + 1] = flatten(node.right)
        return index

    cost_ids = {cost: index for index, cost in enumerate(bucket_indices)}
    roots = array("q", (flatten(index.root) for index in bucket_indices.values()))
    offsets = array("q", [0])
    compatible = array("q")
    for cost in cost_ids:
        compatible.extend(cost_ids[source] for source in dominating_cost_keys[cost])
        offsets.append(len(compatible))
    packed = {
        "points": array("d", (coordinate for precision, _ in values for coordinate in precision)),
        "cost_ids": array("q", (cost_ids[cost] for _, cost in values)),
        "bounds": bounds,
        "nodes": nodes,
        "leaves": leaves,
        "roots": roots,
        "offsets": offsets,
        "compatible": compatible,
    }
    fd = os.memfd_create("approx_ops_prune", os.MFD_CLOEXEC)
    with os.fdopen(fd, "w+b") as shared_file:
        sections = []
        for name, data in packed.items():
            sections.append((name, shared_file.tell(), len(data), data.typecode, data.itemsize))
            data.tofile(shared_file)
        shared_file.flush()
        packed.clear()
        # Workers open this descriptor through /proc; route objects never leave the parent.
        yield (f"/proc/{os.getpid()}/fd/{fd}", tuple(sections))


@contextmanager
def _open_depth_index(descriptor):
    path, sections = descriptor
    with open(path, "rb") as shared_file, mmap.mmap(shared_file.fileno(), 0, access=mmap.ACCESS_READ) as mapping:
        views = {}
        try:
            for name, offset, count, typecode, itemsize in sections:
                views[name] = memoryview(mapping)[offset:offset + count * itemsize].cast(typecode)
            yield views
        finally:
            for view in views.values():
                view.release()


def _shared_has_dominator(index, root, target, exclude_index, require_point_strict):
    bounds, nodes, leaves, points = (index[name] for name in ("bounds", "nodes", "leaves", "points"))
    pending = [root]
    while pending:
        node = pending.pop()
        minimum = tuple(bounds[node * 8:node * 8 + 4])
        if any(value > limit for value, limit in zip(minimum, target)):
            continue
        maximum = tuple(bounds[node * 8 + 4:node * 8 + 8])
        left, right, start, stop, size = nodes[node * 5:node * 5 + 5]
        if (
            not (exclude_index is not None and size == 1)
            and all(value <= limit for value, limit in zip(maximum, target))
            and (not require_point_strict or minimum != target)
        ):
            return True
        if start >= 0:
            for position in range(start, stop):
                point_index = leaves[position]
                if point_index == exclude_index:
                    continue
                point = tuple(points[point_index * 4:point_index * 4 + 4])
                if all(value <= limit for value, limit in zip(point, target)):
                    if not require_point_strict or point != target:
                        return True
        else:
            if tuple(bounds[right * 8:right * 8 + 4]) < tuple(bounds[left * 8:left * 8 + 4]):
                left, right = right, left
            pending.extend((right, left))
    return False


def _dominance_query_batch(descriptor, start, stop):
    dominated = []
    with _open_depth_index(descriptor) as index:
        points, cost_ids, offsets, compatible, roots = (
            index[name] for name in ("points", "cost_ids", "offsets", "compatible", "roots")
        )
        for target_index in range(start, stop):
            target = tuple(points[target_index * 4:target_index * 4 + 4])
            target_cost = cost_ids[target_index]
            for position in range(offsets[target_cost], offsets[target_cost + 1]):
                source_cost = compatible[position]
                same_cost = source_cost == target_cost
                if _shared_has_dominator(
                    index, roots[source_cost], target,
                    target_index if same_cost else None, same_cost,
                ):
                    dominated.append(target_index)
                    break
    return dominated


def _parallel_dominance_mask(executor, descriptor, count, workers, chunk_size):
    chunk_size = chunk_size or max(512, min(16_384, math.ceil(count / (workers * 4))))
    if chunk_size < 1:
        raise ValueError("chunk_size must be positive")
    dominated = bytearray(count)
    futures = []
    try:
        for start in range(0, count, chunk_size):
            futures.append(executor.submit(_dominance_query_batch, descriptor, start, min(start + chunk_size, count)))
        for future in futures:
            for index in future.result():
                dominated[index] = 1
    finally:
        # Keep the shared descriptor alive even if submission or a worker fails.
        wait(futures)
    return dominated


def prune_with_candidates_parallel(
    routes: list[remezData],
    criteria: str,
    prune_range: tuple[str, ...],
    *,
    executor: ProcessPoolExecutor | None = None,
    workers: int = 2,
    chunk_size: int | None = None,
    pruning_type: str = "prune",
) -> tuple[list[remezData], int]:
    """Query the full dominance index in parallel, retaining original route order."""
    if criteria != "depth":
        raise ValueError(f"Unsupported criteria: {criteria}")
    if pruning_type not in ("prune", "lex"):
        raise ValueError("pruning_type must be one of: prune, lex")
    if workers < 1 or (chunk_size is not None and chunk_size < 1):
        raise ValueError("workers and chunk_size must be positive")
    if workers == 1 or len(routes) < 2:
        return _filter_depth_dominated(routes, prune_range, pruning_type)
    values, bucket_indices, dominating_cost_keys = _prepare_depth_index(routes, prune_range, pruning_type)
    with _shared_depth_index(values, bucket_indices, dominating_cost_keys) as descriptor:
        if executor is None:
            with ProcessPoolExecutor(max_workers=workers, mp_context=get_context("spawn")) as pool:
                dominated = _parallel_dominance_mask(pool, descriptor, len(routes), workers, chunk_size)
        else:
            dominated = _parallel_dominance_mask(executor, descriptor, len(routes), workers, chunk_size)
    pruned_routes = [route for index, route in enumerate(routes) if not dominated[index]]
    return pruned_routes, len(routes) - len(pruned_routes)
