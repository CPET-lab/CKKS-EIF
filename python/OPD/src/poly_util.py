
from itertools import product
from math import log2, ceil
from typing import Iterator
from .basic_class import Complexity, Poly, XI, Decomp

def make_all_poly_types(max_deg: int) -> list[list[str]]:
    '''Return every 0/I/F coefficient-type case of exact degree max_deg.'''
    if max_deg < 1:
        raise ValueError("max_deg must be at least 1.")
    return [list(combo) for combo in product(*([['0', 'I', 'F']] * max_deg + [['I', 'F']]))]

def print_dcmp_detail(dcmp: Decomp, step: int):
    print(f"Step {step}")
    print(f"\t"*step, end='')
    print(f"dcmp: {dcmp.restore_dcmp()}")
    print(f"\t"*step, end='')
    print(f"mp: {dcmp.made_powers}")
    print(f"\t"*step, end='')
    dcmp.comp.print_params()
    if dcmp.xi.n != 0:
        print(f"\t"*step, end='')
        dcmp.xi.print_params()
    
    if dcmp.dcmp_p is not None:
        print_dcmp_detail(dcmp.dcmp_p, step+1)
    if dcmp.dcmp_q is not None:
        print_dcmp_detail(dcmp.dcmp_q, step+1)


def solve_xn_routes(multA, n, made_powers=None, minimize_cmult=True) -> list[XI]:
    """
    기능:
        기존에 생성한 x^i를 활용하여 x^n 또는 ax^n을 최소 depth로 구성하는
        경로를 탐색한다. 필요하면 최소 CMult 가지치기를 적용하지 않고 같은
        최소 depth를 만족하는 모든 경로를 반환한다.
    입력:
        multA: False이면 x^n, True이면 ax^n 생성 여부.
        n: 생성할 목표 차수.
        made_powers: 현재 재사용 가능한 순수 x^i 차수 집합.
        minimize_cmult: True이면 최소 CMult 경로만, False이면 최소 depth의
            모든 경로를 유지한다.
    출력:
        선택한 가지치기 조건을 만족하는 XI 객체 목록.

    Rules:
        - made_powers는 순수 power x^k만 추적한다.
        - coefficient route는 순수 x^(a+b)를 만들지 않으므로 made_powers에 추가하지 않는다.
        - multA=True이고 x^n이 이미 존재하면 a * x^n도 고려한다.
        - multA=True이고 x^n이 아직 없더라도, x^n을 만든 뒤 a * x^n 하는 fallback도 고려한다.
        - minimize_cmult=True이면 최소 depth 후보 중 CMult가 최소인 경로만 남긴다.
        - minimize_cmult=False이면 최소 depth를 만족하는 모든 경로를 남긴다.
        - 최종적으로 같은 route/depth/pmult/cmult 후보가 있으면 made_powers가 더 큰 후보만 남긴다.
    """
    from math import ceil, log2

    BASE_COEFF_POWER = 2

    def ceil_log2(v: int) -> int:
        if v <= 1:
            return 0
        return ceil(log2(v))

    def pure_depth(k: int) -> int:
        return ceil_log2(k)

    def direct_coeff_depth(k: int) -> int:

        return ceil_log2(k) + 1

    def base_coeff_depth(k: int) -> int:

        return ceil_log2(k + 1)


    input_made_powers = set(made_powers or set())
    base_made_powers = {0, 1} | input_made_powers


    available_pure = {
        p for p in base_made_powers
        if 1 <= p <= n
    }


    available_coeff_base = set()
    if multA and BASE_COEFF_POWER in base_made_powers and BASE_COEFF_POWER <= n:
        available_coeff_base.add(BASE_COEFF_POWER)


    if not multA and n in base_made_powers:
        xi_obj = XI(multA, n)
        xi_obj.add_routes(
            [],
            pure_depth(n),
            0,
            0,
            set(base_made_powers),
            route_ops=[],
        )
        return [xi_obj]


    min_depths = {}

    for i in range(1, n + 1):

        d_pure = pure_depth(i) if i in available_pure else float("inf")

        for j in range(1, i // 2 + 1):
            k = i - j

            if j in min_depths and k in min_depths:
                p_left = min_depths[j][0]
                p_right = min_depths[k][0]

                if p_left != float("inf") and p_right != float("inf"):
                    p_val = max(p_left, p_right) + 1
                    if p_val < d_pure:
                        d_pure = p_val


        d_coeff = float("inf")

        if multA:

            if i == 1:
                d_coeff = min(d_coeff, 1)


            if i in available_coeff_base:
                d_coeff = min(d_coeff, base_coeff_depth(i))


            if i == n and i in available_pure:
                d_coeff = min(d_coeff, direct_coeff_depth(i))


            if i == n and i not in available_pure and d_pure != float("inf"):
                d_coeff = min(d_coeff, d_pure + 1)

            for j in range(1, i // 2 + 1):
                k = i - j

                if j in min_depths and k in min_depths:

                    c_left = min_depths[j][1]
                    p_right = min_depths[k][0]

                    if c_left != float("inf") and p_right != float("inf"):
                        c_val = max(c_left, p_right) + 1
                        if c_val < d_coeff:
                            d_coeff = c_val


                    if j != k:
                        c_right = min_depths[k][1]
                        p_left = min_depths[j][0]

                        if c_right != float("inf") and p_left != float("inf"):
                            c_val = max(c_right, p_left) + 1
                            if c_val < d_coeff:
                                d_coeff = c_val

        min_depths[i] = (d_pure, d_coeff)


    memo = {}

    def find_paths(target, has_coeff):
        state_key = (target, has_coeff)

        if state_key in memo:
            return memo[state_key]

        target_depth = min_depths[target][1 if has_coeff else 0]

        if target_depth == float("inf"):
            memo[state_key] = []
            return []

        res = []

        if not has_coeff:

            if target in available_pure:
                d = pure_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })

                memo[state_key] = res
                return res


            for j in range(1, target // 2 + 1):
                k = target - j

                if max(min_depths[j][0], min_depths[k][0]) + 1 > target_depth:
                    continue

                l_list = find_paths(j, False)
                r_list = find_paths(k, False)

                for l in l_list:
                    for r in r_list:
                        nd = max(l["depth"], r["depth"]) + 1

                        if nd == target_depth:
                            a, b = sorted((j, k))
                            new_ops = l["ops"] | r["ops"]
                            new_ops.add(("pure", a, b))

                            res.append({
                                "ops": new_ops,
                                "depth": nd,
                            })

        else:

            if target == 1:
                d = 1

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })

                    memo[state_key] = res
                    return res


            if target in available_coeff_base:
                d = base_coeff_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": set(),
                        "depth": d,
                    })


            if target == n and target in available_pure:
                d = direct_coeff_depth(target)

                if d == target_depth:
                    res.append({
                        "ops": {("coeff_direct", 0, target)},
                        "depth": d,
                    })


            if target == n and target not in available_pure:
                pure_target_depth = min_depths[target][0]

                if pure_target_depth != float("inf") and pure_target_depth + 1 == target_depth:
                    pure_paths = find_paths(target, False)

                    for p in pure_paths:
                        if p["depth"] + 1 == target_depth:
                            res.append({
                                "ops": set(p["ops"]),
                                "depth": p["depth"] + 1,
                            })

            for j in range(1, target // 2 + 1):
                k = target - j


                if max(min_depths[j][1], min_depths[k][0]) + 1 <= target_depth:
                    c_list = find_paths(j, True)
                    p_list = find_paths(k, False)

                    for c in c_list:
                        for p in p_list:
                            nd = max(c["depth"], p["depth"]) + 1

                            if nd == target_depth:
                                new_ops = c["ops"] | p["ops"]
                                new_ops.add(("coeff", j, k))

                                res.append({
                                    "ops": new_ops,
                                    "depth": nd,
                                })


                if j != k:
                    if max(min_depths[k][1], min_depths[j][0]) + 1 <= target_depth:
                        c_list = find_paths(k, True)
                        p_list = find_paths(j, False)

                        for c in c_list:
                            for p in p_list:
                                nd = max(c["depth"], p["depth"]) + 1

                                if nd == target_depth:
                                    new_ops = c["ops"] | p["ops"]
                                    new_ops.add(("coeff", k, j))

                                    res.append({
                                        "ops": new_ops,
                                        "depth": nd,
                                    })

        memo[state_key] = res
        return res

    candidates = find_paths(n, multA)


    if not candidates:
        return []

    min_depth = min(c["depth"] for c in candidates)
    candidates = [
        c for c in candidates
        if c["depth"] == min_depth
    ]

    min_ops_count = (
        min(len(c["ops"]) for c in candidates)
        if minimize_cmult
        else None
    )

    kind_order = {
        "pure": 0,
        "coeff_direct": 1,
        "coeff": 2,
    }

    converted_candidates = []

    for cand in candidates:
        ops_set = cand["ops"]

        if min_ops_count is not None and len(ops_set) > min_ops_count:
            continue

        sorted_internal_route = sorted(
            list(ops_set),
            key=lambda op: (
                op[1] + op[2],
                kind_order.get(op[0], 99),
                op[1],
                op[2],
            ),
        )

        updated_made_powers = set(base_made_powers)

        for kind, a, b in sorted_internal_route:

            if kind == "pure":
                updated_made_powers.add(a + b)

        route_ops = [
            (kind, a, b)
            for kind, a, b in sorted_internal_route
        ]
        route = [(a, b) for _, a, b in route_ops]

        depth = cand["depth"]
        pmult = 1 if multA else 0


        cmult = sum(
            1
            for kind, _, _ in sorted_internal_route
            if kind != "coeff_direct"
        )

        converted_candidates.append({
            "route": route,
            "route_ops": route_ops,
            "depth": depth,
            "pmult": pmult,
            "cmult": cmult,
            "made_powers": updated_made_powers,
        })


    pruned_candidates = []

    for idx, item in enumerate(converted_candidates):
        dominated = False

        for jdx, other in enumerate(converted_candidates):
            if idx == jdx:
                continue

            same_external_route = (
                tuple(item["route_ops"]) == tuple(other["route_ops"])
                and item["depth"] == other["depth"]
                and item["pmult"] == other["pmult"]
                and item["cmult"] == other["cmult"]
            )

            if same_external_route and other["made_powers"] > item["made_powers"]:
                dominated = True
                break

        if not dominated:
            pruned_candidates.append(item)

    final_xi_list = []
    seen = set()

    for item in pruned_candidates:
        sig = (
            tuple(item["route_ops"]),
            tuple(sorted(item["made_powers"])),
            item["depth"],
            item["pmult"],
            item["cmult"],
        )

        if sig in seen:
            continue

        seen.add(sig)

        xi_obj = XI(multA, n)
        xi_obj.add_routes(
            item["route"],
            item["depth"],
            item["pmult"],
            item["cmult"],
            item["made_powers"],
            route_ops=item["route_ops"],
        )
        final_xi_list.append(xi_obj)

    return final_xi_list


def _normalize_power_depths(power_cache) -> dict[int, int]:
    '''
    기능:
        set 또는 power-depth dictionary 입력을 실제 depth dictionary로 정규화한다.
    입력:
        power_cache: 순수 power 집합 또는 {power: depth} dictionary.
    출력:
        x^0, x^1을 포함하는 {power: actual depth} dictionary.
    '''
    if isinstance(power_cache, dict):
        result = {
            int(power): int(power_depth)
            for power, power_depth in power_cache.items()
        }
    else:
        result = {
            int(power): ceil(log2(power)) if power > 1 else 0
            for power in set(power_cache or set())
        }
    result.setdefault(0, 0)
    result.setdefault(1, 0)
    return result


def solve_xn_routes_depth_limited(multA: bool, n: int, made_powers: set[int] | dict[int, int] | None = None, max_depth: int | None = None) -> list[XI]:
    '''
    기능:
        입력 power의 실제 depth에서 시작하여 최종 목표 depth가 max_depth 이하인
        x^n 또는 ax^n의 모든 비중복 multiplication DAG를 탐색한다. 중간 power는
        자체 최소 depth를 초과할 수 있다.
    입력:
        multA: False이면 x^n, True이면 ax^n 생성 여부.
        n: 생성할 목표 차수.
        made_powers: 현재 순수 power 집합 또는 각 power의 실제 depth dictionary.
        max_depth: 최종 목표에 허용할 최대 multiplicative depth. None이면 제한하지 않는다.
    출력:
        최종 목표 depth가 max_depth 이하인 서로 다른 XI 객체 목록.
    '''
    from functools import lru_cache

    if n <= 0:
        return []

    base_power_depths = _normalize_power_depths(made_powers)

    @lru_cache(maxsize=None)
    def pure_paths(target: int) -> tuple[frozenset, ...]:
        '''
        기능: 목표 x^target에 기여하는 모든 비중복 pure multiplication DAG를 생성한다.
        입력: target은 생성할 순수 power 차수이다.
        출력: 상세 연산 tuple 집합을 원소로 갖는 tuple이다.
        '''
        paths = set()
        if target in base_power_depths:
            paths.add(frozenset())

        for lhs in range(1, target // 2 + 1):
            rhs = target - lhs
            for left_ops in pure_paths(lhs):
                for right_ops in pure_paths(rhs):
                    paths.add(
                        frozenset(
                            set(left_ops)
                            | set(right_ops)
                            | {("pure", lhs, rhs)}
                        )
                    )
        return tuple(paths)

    @lru_cache(maxsize=None)
    def coeff_paths(target: int) -> tuple[frozenset, ...]:
        '''
        기능: 목표 ax^target에 기여하는 모든 비중복 coefficient DAG를 생성한다.
        입력: target은 생성할 coefficient-bearing power 차수이다.
        출력: 상세 연산 tuple 집합을 원소로 갖는 tuple이다.
        '''
        paths = set()
        for pure_ops in pure_paths(target):
            paths.add(
                frozenset(
                    set(pure_ops)
                    | {("coeff_direct", 0, target)}
                )
            )

        for coeff_power in range(1, target):
            pure_power = target - coeff_power
            for coeff_ops in coeff_paths(coeff_power):
                for pure_ops in pure_paths(pure_power):
                    paths.add(
                        frozenset(
                            set(coeff_ops)
                            | set(pure_ops)
                            | {("coeff", coeff_power, pure_power)}
                        )
                    )
        return tuple(paths)

    def evaluate_ops(ops: frozenset) -> tuple | None:
        '''
        기능: 연산 DAG를 실제 power depth로 평가하고 depth 제한을 검사한다.
        입력: ops는 kind/lhs/rhs 연산으로 구성된 frozenset이다.
        출력: 유효하면 route/depth/cmult/power-depth tuple, 아니면 None이다.
        '''
        kind_order = {
            "pure": 0,
            "coeff_direct": 1,
            "coeff": 2,
        }
        route_ops = sorted(
            ops,
            key=lambda op: (
                op[1] + op[2],
                kind_order.get(op[0], 99),
                op[1],
                op[2],
            ),
        )
        pure_depths = dict(base_power_depths)
        coeff_depths = {}
        generated_pure = set()

        for kind, lhs, rhs in route_ops:
            if kind == "pure":
                target = lhs + rhs
                if lhs not in pure_depths or rhs not in pure_depths:
                    return None
                if target in generated_pure:
                    return None
                generated_depth = max(
                    pure_depths[lhs],
                    pure_depths[rhs],
                ) + 1
                if (
                    target in pure_depths
                    and generated_depth >= pure_depths[target]
                ):
                    return None
                pure_depths[target] = generated_depth
                generated_pure.add(target)
            elif kind == "coeff_direct":
                if rhs not in pure_depths or rhs in coeff_depths:
                    return None
                coeff_depths[rhs] = pure_depths[rhs] + 1
            elif kind == "coeff":
                target = lhs + rhs
                if lhs not in coeff_depths or rhs not in pure_depths:
                    return None
                if target in coeff_depths:
                    return None
                coeff_depths[target] = max(
                    coeff_depths[lhs],
                    pure_depths[rhs],
                ) + 1
            else:
                return None

        target_depths = coeff_depths if multA else pure_depths
        if n not in target_depths:
            return None

        depth = target_depths[n]
        if max_depth is not None and depth > max_depth:
            return None

        cmult = sum(
            1
            for kind, _, _ in route_ops
            if kind != "coeff_direct"
        )
        return route_ops, depth, cmult, pure_depths

    raw_paths = coeff_paths(n) if multA else pure_paths(n)
    results = []
    seen = set()
    for ops in raw_paths:
        evaluated = evaluate_ops(ops)
        if evaluated is None:
            continue
        route_ops, depth, cmult, power_depths = evaluated
        signature = (
            tuple(route_ops),
            depth,
            cmult,
            tuple(sorted(power_depths.items())),
        )
        if signature in seen:
            continue
        seen.add(signature)

        xi = XI(multA, n)
        xi.add_routes(
            [(lhs, rhs) for _, lhs, rhs in route_ops],
            depth,
            1 if multA else 0,
            cmult,
            set(power_depths),
            route_ops=route_ops,
            power_depths=power_depths,
        )
        results.append(xi)

    results.sort(
        key=lambda xi: (
            xi.depth,
            xi.cmult,
            xi.pmult,
            tuple(xi.route_ops),
            tuple(sorted(xi.power_depths.items())),
        )
    )
    return results


def check_without_dcmp(made_powers: set[int], coeff: list[float]) -> bool:
    '''
    Check whether polynomial can be evaluated without cmult.
    '''
    for i, coef in enumerate(coeff):
        if coef != 0:
            if i not in made_powers:
                return False
    return True


def iter_without_dcmp_v2_all(made_powers: set[int] | dict[int, int], coeff: list[float] | list[str], max_depth: int) -> Iterator[tuple[list[tuple[int, int]], list[dict], Complexity, dict[int, int]]]:
    """모든 depth-feasible terminal 계획을 무손실 상태 전이 DAG로 순회한다.

    기존 재귀와 같은 순서로 모든 후보를 생성하지만 각 상태의 Cartesian
    product를 목록에 저장하지 않는다. 상태별로 현재 항의 계획과 다음
    power-depth 상태만 보존하고, 완성된 계획은 DFS 중 하나씩 복원한다.
    """
    poly = Poly(coeff)
    base_power_depths = _normalize_power_depths(made_powers)
    required_powers = [
        degree
        for degree, value in enumerate(poly.coeff)
        if degree >= 1 and value != 0
    ]
    add_count = len(poly.coeff_type) - poly.coeff_type.count("0") - 1
    transitions = {}

    def make_plan(degree: int, coeff_type: str, multA: bool, route_ops) -> dict:
        return {
            "degree": int(degree),
            "coeff_type": str(coeff_type),
            "multA": bool(multA),
            "route_ops": [
                (str(kind), int(lhs), int(rhs))
                for kind, lhs, rhs in route_ops
            ],
        }

    def freeze_plan(plan: dict) -> tuple:
        return (
            plan["degree"],
            plan["coeff_type"],
            plan["multA"],
            tuple(plan["route_ops"]),
        )

    def thaw_plan(plan: tuple) -> dict:
        degree, coeff_type, multA, route_ops = plan
        return make_plan(degree, coeff_type, multA, route_ops)

    def get_transitions(state_key: tuple) -> tuple:
        if state_key in transitions:
            return transitions[state_key]
        pos, power_state = state_key
        current_power_depths = dict(power_state)
        if pos >= len(required_powers):
            transitions[state_key] = ()
            return transitions[state_key]

        degree = required_powers[pos]
        coeff_type = poly.coeff_type[degree]
        term_candidates = []

        if degree in current_power_depths:
            depth = current_power_depths[degree]
            pmult = 1 if coeff_type == "F" else 0
            comp = Complexity()
            comp.insert_value(depth + pmult, 0, pmult, 0)
            if comp.depth <= max_depth:
                term_candidates.append((
                    make_plan(degree, coeff_type, False, []),
                    comp,
                    dict(current_power_depths),
                ))

            if coeff_type == "F":
                for xi in solve_xn_routes_depth_limited(
                    True, degree, current_power_depths, max_depth,
                ):
                    comp = Complexity()
                    comp.insert_value(xi.depth, xi.cmult, xi.pmult, 0)
                    if comp.depth <= max_depth:
                        term_candidates.append((
                            make_plan(degree, coeff_type, True, xi.route_ops),
                            comp,
                            dict(xi.power_depths),
                        ))
        else:
            for xi in solve_xn_routes_depth_limited(
                False, degree, current_power_depths, max_depth,
            ):
                pmult = 1 if coeff_type == "F" else 0
                comp = Complexity()
                comp.insert_value(
                    xi.depth + pmult,
                    xi.cmult,
                    xi.pmult + pmult,
                    0,
                )
                if comp.depth <= max_depth:
                    term_candidates.append((
                        make_plan(degree, coeff_type, False, xi.route_ops),
                        comp,
                        dict(xi.power_depths),
                    ))

            if coeff_type == "F":
                for xi in solve_xn_routes_depth_limited(
                    True, degree, current_power_depths, max_depth,
                ):
                    comp = Complexity()
                    comp.insert_value(xi.depth, xi.cmult, xi.pmult, 0)
                    if comp.depth <= max_depth:
                        term_candidates.append((
                            make_plan(degree, coeff_type, True, xi.route_ops),
                            comp,
                            dict(xi.power_depths),
                        ))

        state_transitions = []
        seen = set()
        for plan, term_comp, term_power_depths in term_candidates:
            next_power_depths = dict(current_power_depths)
            for power, power_depth in term_power_depths.items():
                next_power_depths[power] = min(
                    next_power_depths.get(power, power_depth),
                    power_depth,
                )
            next_state = (
                pos + 1,
                tuple(sorted(next_power_depths.items())),
            )
            transition = (
                freeze_plan(plan),
                (
                    term_comp.depth,
                    term_comp.cmult,
                    term_comp.pmult,
                    term_comp.add,
                ),
                next_state,
            )
            if transition in seen:
                continue
            seen.add(transition)
            state_transitions.append(transition)

        transitions[state_key] = tuple(state_transitions)
        return transitions[state_key]

    initial_state = (0, tuple(sorted(base_power_depths.items())))

    def walk(state_key: tuple, plan_path: tuple, depth: int, cmult: int, pmult: int):
        pos, power_state = state_key
        if pos >= len(required_powers):
            yield plan_path, (depth, cmult, pmult, 0), power_state
            return

        for plan, term_comp, next_state in get_transitions(state_key):
            term_depth, term_cmult, term_pmult, _ = term_comp
            yield from walk(
                next_state,
                plan_path + (plan,),
                max(depth, term_depth),
                cmult + term_cmult,
                pmult + term_pmult,
            )

    for frozen_plans, comp_params, final_power_state in walk(
        initial_state, (), 0, 0, 0,
    ):
        plans = [thaw_plan(plan) for plan in frozen_plans]
        if poly.coeff and poly.coeff[0] != 0:
            plans.insert(
                0,
                make_plan(0, poly.coeff_type[0], False, []),
            )
        comp = Complexity()
        comp.insert_value(*comp_params)
        comp.add = max(add_count, 0)
        routes = [
            (lhs, rhs)
            for plan in plans
            for _, lhs, rhs in plan["route_ops"]
        ]
        yield routes, plans, comp, dict(final_power_state)


def check_without_dcmp_v2_all(made_powers: set[int] | dict[int, int], coeff: list[float] | list[str], max_depth: int) -> list[tuple[list[tuple[int, int]], list[dict], Complexity, dict[int, int]]]:
    """기존 list 반환 API를 유지하는 terminal 완전 탐색 wrapper."""
    return list(iter_without_dcmp_v2_all(made_powers, coeff, max_depth))

def attach(d1: Poly | XI | None, c1: Complexity, d2: Poly, c2: Complexity, attach_type: str) -> Complexity:
    '''
    Attach calculation complexity
    '''
    res = Complexity()
    if d2.coeff == []:
        return c1
    if attach_type == 'x':

        if d2.deg == 0:
            add_mult = 1 if (d2.coeff_type[0] == "F" and d2.coeff[0] != 0) else 0
            res.insert_value(c1.depth + add_mult, c1.cmult, c1.pmult + add_mult, c1.add)

        elif type(d1) == Poly:
            res.depth = max(c1.depth, c2.depth) + 1
            res.cmult = c1.cmult + c2.cmult + 1
            res.pmult = c1.pmult + c2.pmult
            res.add = c1.add + c2.add
        

        elif type(d1) == XI:
            assert d1.n > 0
            res.depth = max(c1.depth, c2.depth) + 1
            res.cmult = c1.cmult + c2.cmult + 1
            res.pmult = c1.pmult + c2.pmult
            res.add = c1.add + c2.add
            
    elif attach_type == '+':
        if d2.deg == 0 and d2.coeff[0] == 0.0:
            return c1
        res.depth = max(c1.depth, c2.depth)
        res.cmult = c1.cmult + c2.cmult
        res.pmult = c1.pmult + c2.pmult
        res.add = c1.add + c2.add + 1
        
    return res
