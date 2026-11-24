# basic.py
from __future__ import annotations
import numpy as np
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from OPS.src.basic_class import remezData, CalData
from OPD.src.algorithm import cal_polyEval
from OPD.src.basic_class import Complexity
import sys
from math import log2, ceil

from common.paths import add_ere_module_path
add_ere_module_path()
from EBModule import EvalStep


# Indicator Function
def evalF_step1(x):
    return 0 if abs(x) >= 0.5 else 1

def evalF_step2(x):
    return 1 if x >= 0.5 else 0

# Interval validation
def verify_interval_step1(interval: tuple) -> bool:
    start1, end1 = interval[0]
    max1 = max(-start1, end1)
    if max1 > interval[1][0]:
        return False
    return True

# Precision calculation
def cal_precision(max_err: float, e_num: int) -> float:
    if max_err <= 0.0:
        return e_num

    if max_err >= 99:
        return -99

    return min(-(log2(max_err)), e_num)

# RouteSearch - 근사 전과 비교하여 중단여부 판단
def check_continue(base: remezData, route: remezData) -> bool:
    # 조건1: 정밀도가 오히려 감소한 경우
    if base.pre_log[-1] > route.pre_log[-1]:
        return False
    # 조건2: 과도한 반복?
    # if len(route.coeff_log) >= 10:
    #     return False
    return True  

# RouteSearch - degree, mode에 따른 더미다항식 생성
def make_dummy_coeff(deg: int, mode: str) -> list[float]:
    return [0.99 if mode == "all" or i % 2 == deg % 2 else 0.0 for i in range(deg + 1)]

# RouteSearch - Remez 근사 전 경로의 예상 실행시간 계산
def precalculate_route_time(
    route: remezData,
    added_depths: set[int],
    time_cache: dict,
) -> dict[int, tuple[float, int]]:
    projected = {}

    for added_depth in added_depths:
        remain_depth = route.total_CalData.depth + added_depth
        prefix_time = 0.0

        # 경로 확장 시 기존 다항식도 증가한 전체 depth에서 다시 평가된다.
        for degree, coeff_type in zip(route.degree_log, route.coeff_type_log):
            prefix_time += time_cache[remain_depth][degree][coeff_type]
            remain_depth -= ceil(log2(degree + 1))

        projected[added_depth] = (prefix_time, remain_depth)

    return projected

# RouteSearch - time cache 검증?
""" def _validate_time_cache_coverage(cache: dict, max_depth: int, max_n: int) -> None:
    for depth in range(2, max_depth + 1):
        if depth not in cache:
            raise ValueError(f"time cache does not contain depth {depth}")

        for degree in range(2, max_n + 1):
            if ceil(log2(degree + 1)) > depth:
                continue

            if degree not in cache[depth]:
                raise ValueError(
                    f"time cache does not contain degree {degree} at depth {depth}"
                )

            expected_types = {
                "all",
                "even" if degree % 2 == 0 else "odd",
            }
            missing_types = expected_types - cache[depth][degree].keys()
            if missing_types:
                raise ValueError(
                    f"time cache is missing {sorted(missing_types)} "
                    f"for degree {degree} at depth {depth}"
                ) 
"""
                    
""" def prune_new_data(orig_data: remezData, new_data: remezData) -> bool:
    '''
    RouteSearch - 정밀도, 연산복잡도 모두 낮은 데이터를 판별
    1. 구간 크기(0/1)
    2. 정밀도(0/1)
    3. 총 연산복잡도
    위 5개의 기준에서 orig_data가 모두 좋다고 판단될 경우 True를 반환함.
    '''
    # 1. 구간 크기 - 작은게 좋다.
    i_orig = orig_data.interval_log[-1]
    i_new = new_data.interval_log[-1]
    orig_size0, orig_size1 = i_orig[0][1] - i_orig[0][0], i_orig[1][1] - i_orig[1][0]
    new_size0, new_size1 = i_new[0][1] - i_new[0][0], i_new[1][1] - i_new[1][0]
    size_flag = bool((orig_size0 < new_size0) and (orig_size1 < new_size1)) 
    
    # 2. 정밀도 - 높은게 좋다
    orig_pre0, orig_pre1 = orig_data.pre_sep_log[-1]
    new_pre0, new_pre1 = new_data.pre_sep_log[-1]
    pre_flag = (orig_pre0 > new_pre0) and (orig_pre1 > new_pre1)
    
    # 3. 총 계산복잡도 - 낮은게 좋다.
    comp_flag = orig_data.total_CalData < new_data.total_CalData

    result = size_flag and pre_flag and comp_flag
    # if size_flag is True and pre_flag is True and comp_flag is True:
    #     print(f"{size_flag}\t{pre_flag}\t{comp_flag}")
    return result """

""" def generate_points(start: np.float64, end: np.float64, step=np.float64(1e-8)):
    if start == end:
        return np.array([start])
    
    max_points = 10000
    while abs(end - start) < step:
        step /= np.float64(10.0)
    step = np.float64(step)

    while True:
        p = np.arange(start, end, step, dtype=np.float64)
        # end 포함 필요 시만 추가
        if p.size == 0 or p[-1] < end:
            p = np.append(p, np.float64(end))
        if p.size <= max_points:
            return p
        step *= np.float64(10.0) """

def generate_points(start: np.float64, end: np.float64, count: int = 10000) -> np.ndarray:
    start = np.float64(start)
    end = np.float64(end)

    if count <= 0:
        raise ValueError("count must be a positive integer")

    if start == end:
        return np.array([start], dtype=np.float64)

    length = abs(end - start)

    # 구간의 길이가 1보다 클 경우 길이 1당 10000개의 점을 샘플링
    if length > 1:
        count = int(np.ceil(length * 10000))

    # 양 끝점을 반드시 포함하기 위해 start != end인 경우 최소 2개 보장
    if count < 2:
        count = 2

    def is_valid(n: int) -> bool:
        pts = np.linspace(start, end, num=n, dtype=np.float64, endpoint=True)
        return np.all(pts[1:] != pts[:-1])

    # 요청한 count개가 가능하면 그대로 사용
    if is_valid(count):
        return np.linspace(start, end, num=count, dtype=np.float64, endpoint=True)

    # 불가능하면 가능한 최대 개수로 줄임
    lo, hi = 2, count
    best = 2

    while lo <= hi:
        mid = (lo + hi) // 2
        if is_valid(mid):
            best = mid
            lo = mid + 1
        else:
            hi = mid - 1

    return np.linspace(start, end, num=best, dtype=np.float64, endpoint=True)

# def generate_points(start: np.float64, end: np.float64, count: int = 10000) -> np.ndarray:
#     start = np.float64(start)
#     end = np.float64(end)

#     if count <= 0:
#         raise ValueError("count must be a positive integer")

#     if start == end:
#         return np.array([start], dtype=np.float64)

#     if count == 1:
#         return np.array([start], dtype=np.float64)

#     def is_valid(n: int) -> bool:
#         pts = np.linspace(start, end, num=n, dtype=np.float64)
#         return np.all(pts[1:] != pts[:-1])

#     # 요청한 count개가 가능하면 그대로 사용
#     if is_valid(count):
#         return np.linspace(start, end, num=count, dtype=np.float64)

#     # 불가능하면 가능한 최대 개수로 줄임
#     lo, hi = 2, count
#     best = 2

#     while lo <= hi:
#         mid = (lo + hi) // 2
#         if is_valid(mid):
#             best = mid
#             lo = mid + 1
#         else:
#             hi = mid - 1

#     return np.linspace(start, end, num=best, dtype=np.float64)

# 함수 종류 확인
def detect_coeff_type(coeff: list[float]) -> str:
    has_even = any(abs(c) != 0.0 for i, c in enumerate(coeff) if i % 2 == 0)
    has_odd  = any(abs(c) != 0.0 for i, c in enumerate(coeff) if i % 2 == 1)

    if has_even and not has_odd:
        return "even"
    if has_odd and not has_even:
        return "odd"
    return "all"

# 특정 함수의 반복횟수에 대한 연산량, 깊이 계산
# def cal_iter(coeff: list[float], iter: int) -> Complexity:
#     res = cal_coeff(coeff)
#     res.depth *= iter
#     res.add_count *= iter
#     res.cmult *= iter
#     res.pmult *= iter

#     return res


# 특정 함수의 연산량, 깊이 계산
# def cal_coeff(coeff: list[float]) -> Complexity:
#     cmp = cal_polyEval(coeff)
#     assert type(cmp) is not bool
#     return cmp.comp
    
#     cal = CalData()
#     cal.depth = cmp.depth
#     cal.cmult = cmp.cmult
#     cal.pmult = cmp.pmult
#     cal.add_count = cmp.add
#     return cal

def cal_complexity(coeff: list[float], es: EvalStep) -> Complexity:
    D, C, P, A = (0, 0, 0, 0)
    if coeff[-1].is_integer():
        D = ceil(log2(len(coeff)-1))
    else:
        D = ceil(log2(len(coeff)))
    
    clist: list[float] = []
    for step in es.eval_step[1:]:
        if step.op == '+':
            A += 1
        else:
            key1 = step.key1[0]
            key1n = int(step.key1[1])
            key2 = step.key2[0]

            if (key1 == 'C' or key2 == 'C'):
                if step.op != 'o':
                    if key1 == 'C' and not clist[key1n].is_integer():
                        P += 1
                    elif key2 == 'C':
                        key2n = int(step.key2[1])
                        if not clist[key2n].is_integer():
                            P += 1 
                else:
                    clist.append(float(step.key2))
                    
            elif (key1 == 'P' and key2 == 'P') or (key1 == 'T' and key2 == 'T'):
                C += 1
            elif (key1 == 'P' and key2 == 'T') or (key1 == 'T' and key2 == 'P'):
                C += 1
            else:
                print(f"cal_complexity ERROR!, key1={step.key1}, key2={step.key2}")
                print("####")
                es.print_step()
                sys.exit()

    comp = Complexity()
    comp.insert_value(D, C, P, A)
    # debug
    # dcmp = cal_polyEval(coeff, {0, 1}, True)
    # comp2 = dcmp.comp
    # if not (comp.depth == comp2.depth and comp.cmult == comp2.cmult and comp.pmult == comp2.pmult and comp.add == comp2.add):
    #     print("COMP ERROR!")
    #     print(coeff)
    #     print(dcmp.restore_dcmp())
    #     comp.print_params()
    #     comp2.print_params()
    #     es.print_step()
    #     sys.exit()
        
    return comp

def cxx_round(x: float) -> int:
    return np.floor(x + 0.5) if x >= 0 else np.ceil(x - 0.5)

def clearNum(a: float) -> str:
    b = cxx_round(a)
    sign = '+' if a > b else '-'

    e = abs(a - float(b))

    if e == 0.0:
        return str(b)

    log2e = float(np.log2(e))

    return f"{b} {sign}2^({log2e:.10f})"

