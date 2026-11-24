from .basic_class import CoefficientNormalizationError, Complexity, Poly, Decomp, XI, NULL_DECOMP
from .poly_util import check_without_dcmp, iter_without_dcmp_v2_all, solve_xn_routes_depth_limited, attach
from .save_class import PowerDepths, PowerState, SearchCache, SearchCacheKey, SearchResult

from math import log2, ceil


"""
cal_polyEval
 ├─ log-depth budget과 입력 power별 실제 depth로 탐색상태 구성
 ├─ 최종 budget을 만족하는 모든 terminal 계획 추가
 ├─ i = 1, ..., deg에 대해
 │   ├─ 중간 power의 최소 depth를 강제하지 않는 모든 x^i 생성 DAG
 │   ├─ 최고차항이 F이면 동일 조건의 모든 ax^i 생성 DAG
 │   └─ f = (x^i 또는 ax^i)p + q로 분해
 │       ├─ p → power cache 확장 → q의 모든 후보 조합
 │       └─ q → power cache 확장 → p의 모든 후보 조합
 ├─ CMult 상한 없이 최종 depth 이하 후보만 유지
 ├─ pruning=True이면 power-depth-aware Pareto frontier만 search cache에 저장
 │  (False이면 구조가 정확히 같은 후보만 제거)
 └─ cal_polyEval wrapper가 CKKS complexity 최소 후보 반환
"""
def cal_polyEval(poly: Poly | list[float] | list[str], made_powers: set[int] | PowerDepths | None = None, is_root: bool = True, divAssumption: bool = True, initial_cache: SearchCache | None = None, intermediate_cache: SearchCache | None = None, depth_budget: int | None = None, cache_root_results: bool = True, pruning: bool = True) -> Decomp:
    """모든 depth-feasible 후보를 순회하고 CKKS complexity 최소 후보를 반환한다.

    ``pruning=True``이면 재귀 후보를 power-depth-aware Pareto frontier로
    축소하고, False이면 기존 전체 후보 반환 방식을 사용한다.
    ``cache_root_results=False``이면 루트 후보 전체를 메모리에 보관하지 않고
    생성 순서대로 최적 후보만 유지한다.
    """
    if cache_root_results:
        candidates = cal_polyEval_candidates(
            poly,
            made_powers,
            depth_budget,
            is_root,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        )
        return candidates[0] if candidates else NULL_DECOMP

    power_depths = _normalize_power_depths(made_powers)
    if not isinstance(poly, Poly):
        poly = Poly(poly)
    if is_root:
        _validate_external_polynomial(poly)
    if depth_budget is None:
        depth_budget = _get_depth_limit(poly)
    if initial_cache is None:
        initial_cache = {}
    if intermediate_cache is None:
        intermediate_cache = {}

    cache = initial_cache if power_depths == {0: 0, 1: 0} else intermediate_cache
    cache_key = _make_cache_key(poly, divAssumption, depth_budget)
    cached_results = _find_cached_decomps(cache, cache_key, power_depths)
    if cached_results:
        return cached_results[0]

    best = None
    for candidate in _iter_generated_candidates(
        poly,
        power_depths,
        depth_budget,
        divAssumption,
        initial_cache,
        intermediate_cache,
        pruning,
    ):
        if best is None or candidate < best:
            best = candidate
    return best if best is not None else NULL_DECOMP


def cal_polyEval_candidates(poly: Poly | list[float] | list[str], made_powers: set[int] | PowerDepths | None = None, depth_budget: int | None = None, is_root: bool = True, divAssumption: bool = True, initial_cache: SearchCache | None = None, intermediate_cache: SearchCache | None = None, pruning: bool = True) -> list[Decomp]:
    '''
    기능:
        중간 power의 개별 최소 깊이를 강제하지 않고 최종 다항식 depth 상한을
        만족하는 모든 비중복 decomposition 및 하위 power 재사용 조합을 반환한다.
    입력:
        poly: 수치 계수, 0/I/F 계수타입 목록 또는 Poly 객체.
        made_powers: 현재 순수 x^i 집합 또는 각 power의 실제 depth dictionary.
        depth_budget: 현재 다항식 평가에 허용되는 최대 multiplicative depth.
        is_root: 초기 다항식 탐색 여부.
        divAssumption: 타입 입력에서 True이면 ax^i 정규화 후 비영 하위 계수를 F로 취급한다.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
        pruning: True이면 안전하게 지배되는 후보를 생성 즉시 제거한다.
    출력:
        pruning이면 Pareto frontier, 아니면 모든 고유 Decomp 후보 목록.
    '''
    # 데이터 초기화
    power_depths = _normalize_power_depths(made_powers)
    made_powers = set(power_depths)

    if not isinstance(poly, Poly):
        poly = Poly(poly)

    if is_root:
        _validate_external_polynomial(poly)
    if depth_budget is None:
        depth_budget = _get_depth_limit(poly)

    # 캐시 설정
    if initial_cache is None:
        initial_cache = {}
    if intermediate_cache is None:
        intermediate_cache = {}

    cache = (
        initial_cache
        if power_depths == {0: 0, 1: 0}
        else intermediate_cache
    )
    cache_key = _make_cache_key(poly, divAssumption, depth_budget)
    
    # key에 해당하는 캐시 후보가 이미 존재하는지 확인
    cached_results = _find_cached_decomps(cache, cache_key, power_depths)
    if len(cached_results) > 0:
        return cached_results

    results = []
    for candidate in _iter_generated_candidates(
        poly,
        power_depths,
        depth_budget,
        divAssumption,
        initial_cache,
        intermediate_cache,
        pruning,
    ):
        _add_decomp_candidate(results, candidate, pruning)

    results = _finish_decomp_candidates(results, pruning)
    results.sort()
    if len(results) > 0:
        _store_cached_decomps(cache, cache_key, power_depths, results)
    return results


def _iter_generated_candidates(poly: Poly, power_depths: PowerDepths, depth_budget: int, divAssumption: bool, initial_cache: SearchCache, intermediate_cache: SearchCache, pruning: bool):
    """논문의 세 OPD case를 기존 순서대로 빠짐없이 생성한다."""
    made_powers = set(power_depths)

    # Case 1: 필요한 power가 모두 cache에 존재한다.
    if check_without_dcmp(made_powers, poly.coeff):
        comp_res = Complexity()
        comp_res.depth = max(
            (power_depths[degree] + (1 if coeff_type == "F" else 0))
            for degree, coeff_type in enumerate(poly.coeff_type)
            if degree >= 1 and coeff_type != "0"
        ) if poly.deg > 0 else 0
        comp_res.cmult = 0
        comp_res.pmult = poly.coeff_type[1:].count("F")
        comp_res.add = max(
            len(poly.coeff_type) - poly.coeff_type.count("0") - 1,
            0,
        )
        if comp_res.depth <= depth_budget:
            term_plans = [
                _make_term_plan(i, ctype, False, [])
                for i, ctype in enumerate(poly.coeff_type)
                if ctype != "0"
            ]
            yield _make_terminal_decomp(
                poly, comp_res, term_plans, power_depths,
            )
    else:
        # Case 2: 분해 없이 필요한 power를 생성한다.
        pinfo = poly.coeff_type if poly.is_type_poly else poly.coeff
        for _, term_plans, comp_res, new_power_depths in iter_without_dcmp_v2_all(
            power_depths, pinfo, depth_budget,
        ):
            yield _make_terminal_decomp(
                poly, comp_res, term_plans, new_power_depths,
            )

    # Case 3: 모든 admissible x^i 및 ax^i 분해를 탐색한다.
    for factor_power in range(1, poly.deg + 1):
        factor_types = [False]
        if poly.coeff_type[poly.deg] == "F":
            factor_types.insert(0, True)

        for multA in factor_types:
            factor_routes = solve_xn_routes_depth_limited(
                multA, factor_power, power_depths, depth_budget,
            )
            for class_xi in factor_routes:
                yield from process_recursion(
                    poly,
                    class_xi,
                    depth_budget,
                    divAssumption,
                    initial_cache,
                    intermediate_cache,
                    pruning,
                )


def process_recursion(poly: Poly, xi: XI, depth_budget: int, divAssumption: bool=True, initial_cache: SearchCache | None = None, intermediate_cache: SearchCache | None = None, pruning: bool = True) -> list[Decomp]:
    '''
    기능:
        하나의 XI에 대한 모든 decomposition 조합을 탐색한다.
    입력:
        poly: 분리 전 부모 다항식.
        xi: 현재 factor의 power 생성정보.
        depth_budget: 부모 다항식에 허용되는 최대 multiplicative depth.
        divAssumption: 타입 입력에서 True이면 ax^i 정규화 후 비영 하위 계수를 F로 취급한다.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
    출력:
        현재 factor에서 생성된 모든 고유 Decomp 후보 목록.
    '''
    try:
        cases = poly.seperate_cases(xi.n, xi.multA, divAssumption)
    except CoefficientNormalizationError:
        return []

    results = []
    for poly_p, poly_q in cases:
        for candidate in _process_recursion_case(
            poly,
            xi,
            poly_p,
            poly_q,
            depth_budget,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        ):
            _add_decomp_candidate(results, candidate, pruning)
    return _finish_decomp_candidates(results, pruning)


def _process_recursion_case(poly: Poly, xi: XI, poly_p: Poly, poly_q: Poly, depth_budget: int, divAssumption: bool=True, initial_cache: SearchCache | None = None, intermediate_cache: SearchCache | None = None, pruning: bool = True) -> list[Decomp]:
    '''
    기능:
        확정된 quotient/remainder 타입에서 p->q 및 q->p의 모든 후보 조합을
        만들고 depth_budget을 만족하는 결과만 반환한다.
    입력:
        poly: 분리 전 부모 다항식.
        xi: 부모와 quotient를 결합할 x^i 또는 ax^i factor.
        poly_p: quotient 다항식.
        poly_q: remainder 다항식.
        depth_budget: 부모 다항식에 허용되는 최대 multiplicative depth.
        divAssumption: 하위 타입 입력의 ax^i 정규화 후 계수타입 처리 방식.
        initial_cache: 입력 power cache가 {0,1}인 후보 cache.
        intermediate_cache: 그 외 입력 power cache의 후보 cache.
    출력:
        depth 조건을 만족하는 모든 고유 Decomp 후보 목록.
    '''
    comp_i = Complexity(xi)
    if comp_i.depth > depth_budget:
        return []

    if poly_p.is_empty() and poly_q.is_empty():
        result = Decomp(poly.coeff, comp_i)
        result.update(xi, None, None)
        return [result]

    if not poly_p.is_empty() and poly_q.is_empty():
        results = []
        for decomp_p in cal_polyEval_candidates(
            poly_p,
            xi.power_depths,
            depth_budget,
            False,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        ):
            comp_total = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, None)
                result.eval_order = "p_only"
                _add_decomp_candidate(results, result, pruning)
        return _finish_decomp_candidates(results, pruning)

    if poly_p.is_empty() and not poly_q.is_empty():
        results = []
        for decomp_q in cal_polyEval_candidates(
            poly_q,
            xi.power_depths,
            depth_budget,
            False,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        ):
            comp_total = attach(xi, comp_i, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, None, decomp_q)
                result.eval_order = "q_only"
                _add_decomp_candidate(results, result, pruning)
        return _finish_decomp_candidates(results, pruning)

    results = []

    for decomp_p in cal_polyEval_candidates(
        poly_p,
        xi.power_depths,
        depth_budget,
        False,
        divAssumption,
        initial_cache,
        intermediate_cache,
        pruning,
    ):
        powers_after_p = _merge_power_depths(
            xi.power_depths,
            decomp_p.merge_power_depths(),
        )
        for decomp_q in cal_polyEval_candidates(
            poly_q,
            powers_after_p,
            depth_budget,
            False,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        ):
            comp_xp = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_total = attach(None, comp_xp, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "p_then_q"
                _add_decomp_candidate(results, result, pruning)

    for decomp_q in cal_polyEval_candidates(
        poly_q,
        xi.power_depths,
        depth_budget,
        False,
        divAssumption,
        initial_cache,
        intermediate_cache,
        pruning,
    ):
        powers_after_q = _merge_power_depths(
            xi.power_depths,
            decomp_q.merge_power_depths(),
        )
        for decomp_p in cal_polyEval_candidates(
            poly_p,
            powers_after_q,
            depth_budget,
            False,
            divAssumption,
            initial_cache,
            intermediate_cache,
            pruning,
        ):
            comp_xp = attach(xi, comp_i, poly_p, decomp_p.comp, 'x')
            comp_total = attach(None, comp_xp, poly_q, decomp_q.comp, '+')
            if comp_total.depth <= depth_budget:
                result = Decomp(poly.coeff, comp_total)
                result.update(xi, decomp_p, decomp_q)
                result.eval_order = "q_then_p"
                _add_decomp_candidate(results, result, pruning)

    return _finish_decomp_candidates(results, pruning)


def _validate_external_polynomial(poly: Poly) -> None:
    if not poly.coeff:
        raise ValueError("At least one coefficient is required.")
    if poly.coeff[-1] == 0:
        raise ValueError("The last coefficient must be nonzero; trailing zeros are not accepted.")


# 다항식 차수와 최고차항 타입으로 depth 상한 계산
def _get_depth_limit(poly: Poly) -> int:
    '''
    기능:
        다항식 차수와 최고차항 타입으로 log-depth 상한을 계산한다.
    입력:
        poly: depth 상한을 계산할 Poly 객체.
    출력:
        최고차항이 F이면 한 단계를 반영한 최대 multiplicative depth.
    '''
    if poly.deg <= 0:
        return 0
    val = 1 if poly.coeff_type[poly.deg] == "F" else 0
    return ceil(log2(poly.deg + val))


def _make_cache_key(poly: Poly, divAssumption: bool, depth_budget: int) -> SearchCacheKey:
    '''
    기능:
        OPD 탐색상태의 다항식 식별자와 탐색옵션을 캐시 키로 변환한다.
    입력:
        poly: 현재 탐색할 Poly 객체.
        divAssumption: 타입 입력의 ax^i 정규화 후 계수타입 처리 방식.
        depth_budget: 현재 탐색에 허용되는 최대 multiplicative depth.
    출력:
        타입 입력은 계수타입, 수치 입력은 실제 계수와 탐색옵션으로 구성된 tuple.
    '''
    poly_key = (
        tuple(poly.coeff_type)
        if poly.is_type_poly
        else tuple(poly.coeff)
    )

    return SearchCacheKey(
        poly_key=poly_key,
        is_type_poly=bool(poly.is_type_poly),
        divAssumption=bool(divAssumption),
        depth_budget=int(depth_budget),
    )


def _find_cached_decomps(
    cache: SearchCache,
    key: SearchCacheKey,
    power_depths: PowerDepths,
) -> list[Decomp]:
    '''
    기능:
        동일한 다항식 타입과 입력 power cache에 저장된 최적 decomposition을 찾는다.
    입력:
        cache: SearchResult tuple 목록을 저장한 dictionary.
        key: 다항식 타입과 탐색옵션으로 구성된 캐시 키.
        power_depths: 현재 입력 pure power의 실제 depth cache.
    출력:
        일치하는 모든 Decomp 후보 목록.
    '''
    power_key = _power_state_key(power_depths)
    candidates = [
        dcmp for powers, dcmp, _ in cache.get(key, [])
        if powers == power_key
    ]

    candidates.sort()
    return candidates


# Store the decomposition candidates in the cache as (power cache, decomposition, complexity) tuples.
def _store_cached_decomps(
    cache: SearchCache,
    key: SearchCacheKey,
    power_depths: PowerDepths,
    results: list[Decomp],
) -> None:
    '''
    기능:
        탐색 후보들을 (power cache, decomposition, complexity) tuple로 저장한다.
    입력:
        cache: SearchResult tuple 목록을 저장할 dictionary.
        key: 다항식 타입과 탐색옵션으로 구성된 캐시 키.
        power_depths: 후보 탐색에 사용된 입력 pure power의 실제 depth cache.
        results: 저장할 Decomp 후보 목록.
    출력:
        없음. cache dictionary를 직접 갱신한다.
    '''
    power_key = _power_state_key(power_depths)
    stored_results = cache.setdefault(key, [])
    stored_signatures = {
        (powers, _make_decomp_signature(dcmp), comp.return_params())
        for powers, dcmp, comp in stored_results
    }

    for result in results:
        signature = (
            power_key,
            _make_decomp_signature(result),
            result.comp.return_params(),
        )

        if signature in stored_signatures:
            continue

        stored_results.append(
            SearchResult(
                power_state=power_key,
                decomposition=result,
                complexity=result.comp,
            )
        )
        stored_signatures.add(signature)


def _power_state_key(power_depths: PowerDepths) -> PowerState:
    '''
    기능:
        power별 실제 depth 상태를 cache 비교용 hashable 값으로 변환한다.
    입력:
        power_depths: {power: actual depth} dictionary.
    출력:
        (power, depth) 쌍으로 구성된 frozenset.
    '''
    return frozenset(
        (int(power), int(power_depth))
        for power, power_depth in power_depths.items()
    )


def _complexity_components(comp: Complexity) -> tuple[int, int, int, int]:
    """Return component-wise costs used by safe Pareto dominance."""
    return comp.depth, comp.cmult, comp.pmult, comp.add


def _power_depths_cover(stronger: PowerDepths, weaker: PowerDepths) -> bool:
    """Whether stronger supplies every weaker power at no greater depth."""
    return all(
        power in stronger and stronger[power] <= power_depth
        for power, power_depth in weaker.items()
    )


def _decomp_dominates(
    candidate: Decomp,
    candidate_powers: PowerDepths,
    other: Decomp,
    other_powers: PowerDepths,
) -> bool:
    """Return True when other cannot beat candidate in any later context."""
    candidate_cost = _complexity_components(candidate.comp)
    other_cost = _complexity_components(other.comp)
    return (
        all(left <= right for left, right in zip(candidate_cost, other_cost))
        and _power_depths_cover(candidate_powers, other_powers)
    )


def _add_decomp_candidate(
    frontier: list[Decomp],
    candidate: Decomp,
    pruning: bool,
) -> None:
    """Append normally or update a power-depth-aware Pareto frontier in place."""
    if not pruning:
        frontier.append(candidate)
        return

    candidate_powers = candidate.merge_power_depths()
    dominated_indices = []
    for index, existing in enumerate(frontier):
        existing_powers = existing.merge_power_depths()
        existing_dominates = _decomp_dominates(
            existing,
            existing_powers,
            candidate,
            candidate_powers,
        )
        candidate_dominates = _decomp_dominates(
            candidate,
            candidate_powers,
            existing,
            existing_powers,
        )

        if existing_dominates:
            # Equivalent resource/cost states need one representative. Match
            # the existing Decomp ordering when it distinguishes the tie.
            if candidate_dominates and candidate < existing:
                frontier[index] = candidate
            return
        if candidate_dominates:
            dominated_indices.append(index)

    for index in reversed(dominated_indices):
        del frontier[index]
    frontier.append(candidate)


def _finish_decomp_candidates(
    results: list[Decomp],
    pruning: bool,
) -> list[Decomp]:
    """Preserve the legacy exact-structure dedupe when pruning is disabled."""
    return results if pruning else _dedupe_decomps(results)


# 중복된 후보 제거
def _dedupe_decomps(results: list[Decomp]) -> list[Decomp]:
    '''
    기능:
        구조와 complexity가 완전히 같은 decomposition 후보만 제거한다.
    입력:
        results: 중복을 포함할 수 있는 Decomp 후보 목록.
    출력:
        최초 등장 순서를 유지한 고유 Decomp 후보 목록.
    '''
    unique_results = []
    seen = set()
    for result in results:
        signature = (
            _make_decomp_signature(result),
            result.comp.return_params(),
            tuple(sorted(result.merge_power_depths().items())),
        )
        if signature in seen:
            continue
        seen.add(signature)
        unique_results.append(result)
    return unique_results


def _make_decomp_signature(dcmp: Decomp | None) -> tuple | None:
    '''
    기능:
        중복 캐시 후보를 판별할 수 있도록 decomposition 트리를 tuple로 변환한다.
    입력:
        dcmp: signature를 생성할 Decomp 객체 또는 None.
    출력:
        계수타입, factor, 평가순서, 자식구조를 포함한 tuple 또는 None.
    '''
    if dcmp is None:
        return None

    return (
        tuple(dcmp.coeff_type),
        bool(dcmp.xi.multA),
        int(dcmp.xi.n),
        tuple(dcmp.xi.route_ops),
        tuple(sorted(dcmp.xi.power_depths.items())),
        str(dcmp.eval_order),
        tuple(
            (
                int(plan["degree"]),
                str(plan["coeff_type"]),
                bool(plan["multA"]),
                tuple(plan["route_ops"]),
            )
            for plan in dcmp.term_plans
        ),
        _make_decomp_signature(dcmp.dcmp_p),
        _make_decomp_signature(dcmp.dcmp_q),
    )


def _make_terminal_decomp(
    block: Poly,
    comp: Complexity,
    term_plans: list[dict],
    power_depths: set[int] | dict[int, int],
) -> Decomp:
    '''
    기능:
        항별 terminal 평가계획과 출력 power-depth cache를 leaf Decomp로 구성한다.
    입력:
        block: terminal로 평가한 다항식.
        comp: terminal 평가의 계산복잡도.
        term_plans: 항별 power 및 계수 결합 계획.
        power_depths: 평가 후 사용 가능한 순수 power의 실제 depth cache.
    출력:
        terminal 평가계획을 저장한 Decomp 객체.
    '''
    power_depths = _normalize_power_depths(power_depths)
    route_ops = [
        op
        for plan in term_plans
        for op in plan["route_ops"]
    ]
    routes = [(a, b) for _, a, b in route_ops]

    xi = XI(False, 0)
    xi.add_routes(
        routes,
        comp.depth,
        comp.pmult,
        comp.cmult,
        set(power_depths),
        route_ops=route_ops,
        power_depths=power_depths,
    )

    result = Decomp(block.coeff, comp, xi)
    result.eval_order = "terminal"
    result.term_plans = [
        {
            "degree": plan["degree"],
            "coeff_type": plan["coeff_type"],
            "multA": plan["multA"],
            "route_ops": list(plan["route_ops"]),
        }
        for plan in term_plans
    ]
    return result


def _append_constant_term_plan(poly: Poly, term_plans: list[dict]) -> None:
    if poly.coeff and poly.coeff[0] != 0:
        term_plans.insert(
            0,
            _make_term_plan(
                degree=0,
                coeff_type=poly.coeff_type[0],
                multA=False,
                route_ops=[],
            ),
        )


def _make_term_plan(degree: int, coeff_type: str, multA: bool, route_ops) -> dict:
    return {
        "degree": int(degree),
        "coeff_type": str(coeff_type),
        "multA": bool(multA),
        "route_ops": _copy_route_ops(route_ops),
    }


def _copy_route_ops(route_ops) -> list[tuple[str, int, int]]:
    return [
        (str(kind), int(a), int(b))
        for kind, a, b in route_ops
    ]


def _normalize_power_depths(power_cache) -> PowerDepths:
    '''
    기능:
        set 또는 power-depth dictionary를 OPD가 사용하는 실제 depth 상태로 변환한다.
    입력:
        power_cache: 순수 power 집합, {power: depth} dictionary 또는 None.
    출력:
        x^0, x^1을 포함하는 power-depth dictionary.
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


def _merge_power_depths(*power_caches: PowerDepths) -> PowerDepths:
    '''
    기능:
        여러 power-depth cache를 병합하고 같은 power에는 가장 작은 depth를 유지한다.
    입력:
        power_caches: 병합할 power-depth dictionary 목록.
    출력:
        병합된 power-depth dictionary.
    '''
    result = {0: 0, 1: 0}
    for power_cache in power_caches:
        for power, power_depth in power_cache.items():
            result[power] = min(
                result.get(power, power_depth),
                power_depth,
            )
    return result


def _ceil_log2(v: int) -> int:
    if v <= 1:
        return 0
    return ceil(log2(v))
