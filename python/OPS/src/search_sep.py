# Remez route search algorithm(BFS)
from math import ceil, log2
import copy, sys, time
from tqdm import tqdm
from types import SimpleNamespace
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import get_context
from contextlib import redirect_stdout, redirect_stderr, ExitStack

from common.remez.algorithm import remez_algorithm
from common.remez.print import print_intervals, coeff2txt

from OPS.src.basic_class import remezData, CalData
from OPS.src.basic import *
from OPS.src.write_progress import TqdmOverwriteFile, DiskRouteSet, RouteChunkWriter, _save_routes_to_disk, _iter_disk_routes, cleanup_stale_remez_tempdirs, create_remez_tempdir
from common.runtime_cache import parse_time_cache
from OPS.src.prune import PARALLEL_MIN_ROUTES, prune_with_BEST, prune_with_candidates

# Boundary Calculation
from common.paths import OPD_CACHE, RUNTIME_CACHE, OPS_COEFFICIENTS, add_ere_module_path
add_ere_module_path()
import EBModule as ebm
from EBModule import ErrBound, EvalStep, estimate_range

import os, shutil

ROUTE_SPILL_THRESHOLD = 1_000_000
ROUTE_CHUNK_SIZE = 100_000
FUTURES_PER_WORKER = 4

PRUNE_RANGE_ORDER = ("c", "p", "a")

CACHEFILE = str(OPD_CACHE)
WORKER_CACHE = None

time_cache = None
eb = None
key_func = None

def _init_worker(cachefile: str, eb_params: tuple, time_cache_file: str):
    global WORKER_CACHE, eb, time_cache
    WORKER_CACHE = ebm.DecompCache(cachefile)
    eb = ErrBound(*eb_params)
    time_cache = parse_time_cache(time_cache_file)

def _get_worker_cache():
    global WORKER_CACHE
    if WORKER_CACHE is None:
        raise RuntimeError("WORKER_CACHE is not initialized")
    return WORKER_CACHE

def multi_remez(options: SimpleNamespace) -> None:
    global time_cache, eb, key_func
    eb = ErrBound(*options.eb_params)

    if options.criteria == "depth" and options.pruning_type not in ("base", "prune", "lex"):
        raise ValueError("pruning_type must be one of: base, prune, lex")

    prune_range = tuple(getattr(options, "prune_range", ()))
    if prune_range != PRUNE_RANGE_ORDER[:len(prune_range)]:
        raise ValueError(
            "prune_range must be one of [], ['c'], ['c', 'p'], "
            "or ['c', 'p', 'a']."
        )
    options.prune_range = prune_range
    
    if options.criteria == "depth":
        key_func = lambda route: (
            route.total_CalData.key("depth"),
            tuple(cd.key("depth") for cd in route.CalData_log),
        )
    elif options.criteria == "time":
        key_func = lambda route: route.total_time
    
    time_cache = parse_time_cache(RUNTIME_CACHE)
    options.max_search_depth = getattr(options, "max_search_depth", max(time_cache))
    if options.max_search_depth > max(time_cache):
        raise ValueError(
            f"max_search_depth={options.max_search_depth} exceeds the time cache "
            f"maximum depth {max(time_cache)}"
        )
    
    with ExitStack() as stack:
        if options.FILENAME:
            out_file = stack.enter_context(open(options.FILENAME, "w", encoding="utf-8"))
            stack.enter_context(redirect_stdout(out_file))
            stack.enter_context(redirect_stderr(out_file))
        
        if options.criteria == "depth" and options.pruning_type in ("prune", "lex"):
            attributes = ("precision0(lb,ub)", "precision1(lb,ub)", "d", *options.prune_range)
            print(f"Prune attributes: {', '.join(attributes)}")
            if options.pruning_type == "lex":
                print("Prune comparison: lexicographic")

        cleanup_stale_remez_tempdirs() # 임시파일 제거
        
        for p_num in range(options.start, options.end+1):
            options.BEST = None
            p = pow(2, p_num)
            e = eb.Bc / eb.scale
        
            # 초기 구간 설정
            if (p_num == 5 and options.e_num == 15) or (p_num == 8 and options.e_num == 30):
                options.init_intervals = ((0, e), (1-e, p))
            else:
                options.init_intervals = ((0, e), (1-e, p-1+e))
            print(f"p=2^{p_num}, e=2^{options.e_num}")
        
            # Step1 - 각 n에 대하여 초기 근사를 진행, 각 depth에 대해 기본 객체 생성
            routes_step1: list[remezData] = remez_step1(options)

            # Step2 - criteria에 대하여 최적의 결과 산출
            BEST_ROUTE: remezData | None = remez_step2(routes_step1, options)
            if not BEST_ROUTE:
                print("Approximation Failed: No valid route generated.")
                print("\n##########\n")
                continue

            # Step3 - 저장
            print(f"(p, e): {p_num}, {options.e_num}")
            print(f"Min {options.criteria}")
            BEST_ROUTE.print_params()
            if options.save_result:
                filename = OPS_COEFFICIENTS / f"coeff_{p_num}_{options.e_num}_{options.criteria}.txt"
                filename.parent.mkdir(parents=True, exist_ok=True)
                with open(filename, "w", encoding="utf-8") as f:
                    for coeff in BEST_ROUTE.coeff_log:
                        f.write(coeff2txt(coeff))
                        f.write('\n')
                
            print("\n##########\n")

# Step 1 - Approximate Even Functions
def remez_step1(options: SimpleNamespace) -> list[remezData]:
    routes = []
    cache = ebm.DecompCache(CACHEFILE)
    
    max_n = options.max_n
    intervals = options.init_intervals
    
    for n in range(max_n - (max_n % 2), 1, -2):
        try:
            coeff = remez_algorithm(n, intervals, 1, "even", "normal")
        except Exception as e:
            print(f"Remez error!\n{e}")
            sys.exit()

        if coeff == [-1]:
            print(f"Init approximation\nn={n} failed.")
            continue
        else:
            # Check interval disjoint
            es = cache.restore_EvalStep(coeff)
            max_err, pre_sep, next_intervals = calculate_next_remez(es, 1, intervals, options)
            if verify_interval_step1(next_intervals):
                print(f"Init approximation\n n={n}, step 0, intervals: ", end='')
                print_intervals(next_intervals, 10)
            else:
                print(f"Init approximation\nn={n} failed.")
                continue
        print('-' * 30)

        pre = -(log2(max_err))
        data = remezData(intervals)
        data.update(coeff, pre, pre_sep, es, next_intervals, "even")
        routes.append(data)
        
    return routes

# Step 2 - BFS Search
def remez_step2(routes: list[remezData], options: SimpleNamespace):
    global key_func
    step_count = 1
    worker_count = options.cpu_count if options.multiThread else 1
    spill_enabled = options.criteria == "time" or options.pruning_type == "base"
    spill_root = create_remez_tempdir() if spill_enabled else None

    options.BEST = None

    t0 = time.perf_counter()
    try:
        time_cache_file = str(RUNTIME_CACHE)
        with ProcessPoolExecutor(max_workers=worker_count, initializer=_init_worker, initargs=(CACHEFILE, options.eb_params, time_cache_file), mp_context=get_context("spawn")) as executor:
            while len(routes) > 0:
                prev_disk_routes = routes if isinstance(routes, DiskRouteSet) else None

                routes, temp_BEST = step_with_MultiThread(executor, routes, step_count, worker_count, spill_root, options)

                if temp_BEST is not None:
                    if options.BEST is None:
                        options.BEST = temp_BEST
                    else:
                        options.BEST = min([options.BEST, temp_BEST], key=key_func)

                if prev_disk_routes is not None:
                    prev_disk_routes.cleanup()

                step_count += 1

    except Exception as e:
        print(f"Step2 error!\n{e}")
        return None

    finally:
        dt = time.perf_counter() - t0
        print(f"[Step2] total_time={dt:.2f}s")
        if spill_root is not None:
            shutil.rmtree(spill_root, ignore_errors=True)

    return options.BEST

# Search with MultiThreading
def step_with_MultiThread(executor: ProcessPoolExecutor, routes: list[remezData] | DiskRouteSet, step_count: int, worker_count: int, spill_root: str | None, options: SimpleNamespace):
    global key_func

    route_count = len(routes)
    print(f"Step {step_count}: Processing {route_count} routes...")
    step_t0 = time.perf_counter()

    BEST_ROUTE = options.BEST

    total_count = 0
    completed_count = 0
    discard_count = 0
    prune_count1 = 0
    prune_count2 = 0

    if sys.stdout.seekable():
        pbar_file = TqdmOverwriteFile(sys.stdout)
    else:
        pbar_file = sys.stdout

    spill_enabled = options.criteria == "time" or options.pruning_type == "base"
    if spill_enabled and spill_root is None:
        raise RuntimeError("spill_root is required for spill-enabled search")

    step_dir = None
    in_dir = None
    out_dir = None
    input_spilled = False

    if spill_enabled:
        step_dir = os.path.join(spill_root, f"step_{step_count:04d}")
        in_dir = os.path.join(step_dir, "input")
        out_dir = os.path.join(step_dir, "output")

    if spill_enabled and route_count >= ROUTE_SPILL_THRESHOLD:
        os.makedirs(step_dir, exist_ok=True)

        if isinstance(routes, DiskRouteSet):
            routes_disk = routes
        else:
            routes_disk = _save_routes_to_disk(routes, in_dir, ROUTE_CHUNK_SIZE)

        route_iter = _iter_disk_routes(routes_disk)
        input_spilled = True
    else:
        if isinstance(routes, DiskRouteSet):
            loaded_routes = []
            for _, chunk_routes in _iter_disk_routes(routes):
                loaded_routes.extend(chunk_routes)
            routes = loaded_routes

        route_iter = [(None, routes)]

    candidate_buffer: list[remezData] = []
    writer = None
    step_best = None

    futures_chunk_size = max(worker_count * FUTURES_PER_WORKER, 1)

    with tqdm(total=route_count, desc=f"Step {step_count}", file=pbar_file, mininterval=1.0, ascii=True) as pbar:
        for _, chunk_routes in route_iter:
            for i in range(0, len(chunk_routes), futures_chunk_size):
                sub_routes = chunk_routes[i:i + futures_chunk_size]
                futures = [executor.submit(expand_route, route, BEST_ROUTE, options) for route in sub_routes]

                for future in as_completed(futures):
                    local_completed, local_new, calculated_count, discarded_count = future.result()

                    total_count += calculated_count
                    discard_count += discarded_count

                    if local_completed:
                        completed_count += len(local_completed)
                        local_best = min(local_completed, key=key_func)

                        if step_best is None:
                            step_best = local_best
                        else:
                            step_best = min([step_best, local_best], key=key_func)

                    if local_new:
                        start_writer = spill_enabled and writer is None and (input_spilled or len(candidate_buffer) + len(local_new) >= ROUTE_SPILL_THRESHOLD)

                        if start_writer:
                            os.makedirs(out_dir, exist_ok=True)
                            writer = RouteChunkWriter(out_dir, ROUTE_CHUNK_SIZE)

                            for data in candidate_buffer:
                                writer.add(data)

                            candidate_buffer.clear()

                        if writer is not None:
                            for data in local_new:
                                writer.add(data)
                        else:
                            candidate_buffer.extend(local_new)

                    pbar.update(1)

            del chunk_routes

    if writer is not None:
        total_candidate = writer.close()
    else:
        total_candidate = candidate_buffer

    # 완료 데이터 처리
    if step_best is not None:
        if BEST_ROUTE is None:
            BEST_ROUTE = step_best
        else:
            BEST_ROUTE = min([BEST_ROUTE, step_best], key=key_func)

    search_dt = time.perf_counter() - step_t0

    # Pruning
    pruning_t0 = time.perf_counter()
    best_pruning_dt = 0.0
    candidate_pruning_dt = 0.0
    candidate_workers = 1
    if options.criteria == "time" or options.pruning_type in ("base", "prune", "lex"):
        if BEST_ROUTE is not None:
            best_pruning_t0 = time.perf_counter()
            total_candidate, pc1 = prune_with_BEST(total_candidate, BEST_ROUTE, options)
            best_pruning_dt = time.perf_counter() - best_pruning_t0
            prune_count1 += pc1

        assert len(total_candidate) == total_count - completed_count - discard_count - prune_count1, "Route count mismatch after BEST pruning."

    if options.criteria == "depth" and options.pruning_type in ("prune", "lex"):
        candidate_pruning_t0 = time.perf_counter()
        parallel_min_routes = getattr(options, "pruning_min_routes", PARALLEL_MIN_ROUTES)
        if worker_count > 1 and len(total_candidate) >= max(2, parallel_min_routes):
            candidate_workers = worker_count
        total_candidate, pc2 = prune_with_candidates(
            total_candidate,
            options.criteria,
            options.prune_range,
            executor=executor,
            workers=worker_count,
            parallel_min_routes=parallel_min_routes,
            pruning_type=options.pruning_type,
        )
        candidate_pruning_dt = time.perf_counter() - candidate_pruning_t0
        prune_count2 += pc2

    assert len(total_candidate) == total_count - completed_count - discard_count - prune_count1 - prune_count2, "Route count mismatch after candidate pruning."
    pruning_dt = time.perf_counter() - pruning_t0

    print(f"\tTime: search={search_dt:.2f}s, pruning={pruning_dt:.2f}s")
    print(f"\tPruning detail: BEST={best_pruning_dt:.2f}s, candidates={candidate_pruning_dt:.2f}s, candidate_workers={candidate_workers}")
    print(f"\tRoutes: newly calculated={total_count}, completed={completed_count}")

    if options.criteria == "time" or options.pruning_type in ("base", "prune", "lex"):
        print(f"\tPrune(vs. BEST): {prune_count1} routes pruned.")
    if options.criteria == "depth" and options.pruning_type in ("prune", "lex"):
        print(f"\tPrune(vs. cand): {prune_count2} routes pruned.")

    print('-' * 30)

    return total_candidate, BEST_ROUTE

# Step 2 경로 확장 프로세스
def expand_route(route: remezData, BEST: remezData | None, options: SimpleNamespace) -> tuple[list[remezData], list[remezData], int, int]:
    global eb, time_cache
    cache = _get_worker_cache()

    N_LIMIT = options.max_n
    
    local_completed: list[remezData] = []
    local_new: list[remezData] = []
    calculated_count = 0
    discarded_count = 0

    remain_depth = options.max_search_depth - route.total_CalData.depth
    if remain_depth < 2:
        return local_completed, local_new, calculated_count, discarded_count
    N_LIMIT = min(options.max_n, (1 << remain_depth) - 1)
        
    # depth 기준 경로확장 제한
    if options.criteria == "depth" and BEST is not None:
        remain_depth = min(remain_depth, BEST.total_CalData.depth - route.total_CalData.depth)
        if remain_depth < 2:
            return local_completed, local_new, calculated_count, discarded_count
        N_LIMIT = min(options.max_n, (1 << remain_depth) - 1)

    # time 기준에서는 BEST가 있을 때만 예상시간을 상한으로 사용할 수 있다.
    # 실제 Remez 근사 전에 기존 경로가 추가 depth에서 소요할 시간을 계산한다.
    time_precheck = options.criteria == "time" and BEST is not None
    precalculated_time = {}
    if time_precheck:
        added_depths = {ceil(log2(n + 1)) for n in range(2, N_LIMIT + 1)}
        precalculated_time = precalculate_route_time(route, added_depths, time_cache)

    # 경로 확장
    for n in range(2, N_LIMIT + 1):
        modes_to_try = ("even" if n % 2 == 0 else "odd", "all")
        added_depth = ceil(log2(n + 1))
        for mode in modes_to_try:
            try:
                intervals = route.interval_log[-1]

                # 예상시간이 이미 BEST보다 크면 근사 후에도 개선될 수 없으므로
                # 비용이 큰 Remez 근사를 수행하지 않는다.
                if time_precheck:
                    prefix_time, polynomial_depth = precalculated_time[added_depth]
                    projected_time = prefix_time + time_cache[polynomial_depth][n][mode]
                    if projected_time > BEST.total_time:
                        continue

                coeff = remez_algorithm(n, intervals, 2, mode, "normal")

                if coeff == [-1]:
                    continue
                
                # 새로운 치역 계산
                es = cache.restore_EvalStep(coeff)
                max_err, pre_sep, next_intervals = calculate_next_remez(es, 2, intervals, options)

                if max_err == 99:
                    calculated_count += 1
                    discarded_count += 1
                    continue

                pre = min(-(log2(max_err)), options.e_num)

                expanded_data = copy.deepcopy(route)
                expanded_data.update(coeff, pre, pre_sep, es, next_intervals, mode)
                expanded_data.update_time(time_cache)
                calculated_count += 1

                if pre >= options.e_num:
                    local_completed.append(expanded_data)

                elif expanded_data.total_CalData.depth + 2 > options.max_search_depth:
                    discarded_count += 1

                elif route.pre_log[-1] < expanded_data.pre_log[-1]:
                    local_new.append(expanded_data)

                else:
                    discarded_count += 1

            except Exception:
                continue

    # Worker에 전달된 BEST를 기준으로 기존 즉시 제거 조건을 적용한다.
    if (
        BEST is not None
        and local_new
        and (options.criteria == "time" or options.pruning_type in ("base", "prune", "lex"))
    ):
        local_new, pruned_count = prune_with_BEST(local_new, BEST, options, depth_margin=1)
        discarded_count += pruned_count

    assert calculated_count == len(local_completed) + len(local_new) + discarded_count, "Count mismatch in remez_step."

    return local_completed, local_new, calculated_count, discarded_count

# 치역 계산
def calculate_next_remez(es: EvalStep, evalF_mode: int, intervals: tuple, options: SimpleNamespace):
    '''
    다음 근사를 위해 최대오차(max_err), 치역(next_interval) 계산 \\
        @return
        max_err: float: 모든 구간 통합 최대오차 \\
        err0, err1: float: 각 목표값에 대한 정밀도 \\
        next_interval: 치역구간
    '''
    global eb
    next_interval = ()
    max_err = 0.0
    err0 = 0.0
    err1 = 0.0
    assert evalF_mode in [1, 2], "Invalid evalF_mode"
    interval_size = options.interval_size if options.interval_size else pow(2, -14)
    
    if evalF_mode == 1:
        bc = eb.Bc / eb.scale
        next_interval0 = estimate_range(eb, (0.0, bc), es, interval_size)

        p = round(intervals[1][1])
        next_interval1_candidates = []
        for i in range(1, p+1):
            ni = estimate_range(eb, (i - bc, i + bc), es, interval_size)
            next_interval1_candidates.append(ni)
        next_interval1 = (min(x[0] for x in next_interval1_candidates),
                          max(x[1] for x in next_interval1_candidates))
        next_interval = (next_interval1, next_interval0)
    else:
        next_interval0 = estimate_range(eb, intervals[0], es, interval_size)
        next_interval1 = estimate_range(eb, intervals[1], es, interval_size)
        next_interval = (next_interval0, next_interval1)

    for start, end in next_interval:
        rounded_start = round(start)
        rounded_end = round(end)

        if rounded_start != rounded_end or rounded_start not in (0, 1):
            return 99, (99, 99), ((-99, 99), (-99, 99))

        if rounded_start == 0:
            err0 = min(-log2(abs(start)), -log2(abs(end)), options.e_num)
            max_err = max(abs(start), abs(end), max_err)

        else:
            err1 = min(-log2(abs(1-start)), -log2(abs(1-end)), options.e_num)
            max_err = max(abs(1-start), abs(1-end), max_err)
    
    return max_err, (err0, err1), next_interval
