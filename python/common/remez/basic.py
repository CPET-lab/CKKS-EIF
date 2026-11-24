import numpy as np
from numpy.polynomial import Polynomial
from .print import debug_print

''' 기타 필요한 함수들 '''
# 다항함수 평가
def evalP(coeff, x):
    return sum(coeff[i] * x**i for i in range(len(coeff)))

def error_abs(coeff, x, y):
    return abs(evalP(coeff, x) - y)

def error_func(coeff, x, y):
    return evalP(coeff, x) - y

# 구간 [a, b]에서 chebyshev node를 사용하여 m개의 점 샘플링
def sample_points(a: np.float64, b: np.float64, m) -> np.ndarray:
    if m <= 0:
        print("Sample 0 points!!!")
        return np.ndarray([], dtype=float)
    k = np.arange(m)
    x = np.cos(np.pi * (2 * k + 1) / (2 * m))
    return 0.5 * (b + a) + 0.5 * (b - a) * x

def sample_points_multi(m, intervals, print_mode) -> np.ndarray:
    # 구간의 총 길이를 계산
    lengths = np.array([end - start for start, end in intervals])
    total_length = lengths.sum()
    n = len(intervals)
    m_alloc = np.ones(n, dtype=int)
    left = m - n

    if total_length == 0:
        # 총 길이가 0이면 m을 균등하게 분배
        base = m // n
        remainder = m % n
        m_alloc = np.full(n, base, dtype=int)
        if remainder > 0:
            m_alloc[:remainder] += 1
    else:
        # 각 구간의 길이 비율로 m을 분배, 남은 점은 긴 길이부터 추가로 할당
        m_float = left * lengths / total_length
        additional = np.floor(m_float).astype(int)
        m_alloc += additional

        remainder = m - m_alloc.sum()
        if remainder > 0:
            idx_desc = np.argsort(-lengths)
            for idx in idx_desc[:remainder]:
                m_alloc[idx] += 1

    # 구간 별 샘플링
    chunks = []
    for (s, e), num in zip(intervals, m_alloc):
        debug_print(f"sample {num} points in interval[{s}, {e}]", print_mode)
        if num > 0:
            chunks.append(sample_points(s, e, num))

    x = np.concatenate(chunks)
    x.sort()
    return x

# 구간 변경
def slice_interval(approx_mode: str, intervals: tuple) -> tuple:
    if approx_mode == "all":
        return intervals
    
    # 기존 구간 intervals 정렬
    new_intervals = []
    for start, end in intervals:
        if start < 0 and end < 0:
            continue
        elif start < 0 and end > 0:
            if abs(start) > abs(end):
                new_intervals.append((0.0, abs(start)))
            else:
                new_intervals.append((0.0, end))
        else:
            new_intervals.append((start, end))

    return tuple(new_intervals)
            
# 행렬식 구성
def create_matrix(powers: list, x: np.ndarray, evalF) -> tuple[np.ndarray, np.ndarray]:
    A_matrix = []
    y_matrix = []
    for xindex in range(len(x)):
        x_sample = x[xindex]
        rowA = [x_sample ** p for p in powers] + [(-1) ** xindex]
        y = evalF(x_sample)
        A_matrix.append(rowA)
        y_matrix.append(y)

    A_matrix = np.array(A_matrix, dtype=np.float64)
    y_matrix = np.array(y_matrix, dtype=np.float64)
    return A_matrix, y_matrix

# 행렬식 계산 및 계수벡터 연산
def solve_matrix(A: np.ndarray, y: np.ndarray, n: int, powers: list) -> tuple[list, float]:
    B = np.linalg.solve(A, y)
    E = B[-1]
    coeff = []
    for k in range(n + 1):
        if k in powers:
            coeff.append(B[0])
            B = B[1:]
        else:
            coeff.append(0.0)
    return coeff, E


# 지역 극점의 x좌표 계산
def calculate_local_max(coeff, evalF, intervals: tuple) -> tuple:
    poly = Polynomial(coeff)  # 다항식 객체 생성
    deriv_roots = poly.deriv().roots() # 도함수의 근 계산
    real_roots = deriv_roots[np.isreal(deriv_roots)].real # 실근만 추리기
    
    max_point_x = []
    max_point_y = []
    
    # 각 구간에 속하는 점만 선별
    for start, end in intervals:
        valid_roots = real_roots[(real_roots > start) & (real_roots < end)]
        
        # 구간 시작/끝점, 구간 내 극점들
        candidates = np.unique(np.concatenate(([start, end], valid_roots)))
        
        f_vals = np.array([evalF(x) for x in candidates])
        p_vals = poly(candidates)
        errors = p_vals - f_vals
        
        # 극점 모두 추가
        max_point_x.extend(candidates)
        max_point_y.extend(errors)

    return max_point_x, max_point_y

# mu(x) 계산
def eval_mu(coeff, x: float, eps: float = 0.0) -> int:
        poly = Polynomial(coeff)
        r2 = poly.deriv(2)(x)  # p''(x)
        if r2 < -eps:
            return 1
        if r2 > eps:
            return -1
        return 0 
    
# 논문에서 제시한 극점 선별 알고리즘 구현(Algorithm 3)
def select_max_points(max_point_x: list, max_point_y: list, sample_number: int, coeff: list, evalF, print_mode: str = "normal", E: float | None = None):
    errs = []
    for x, y in zip(max_point_x, max_point_y):
        err = evalP(coeff, x) - evalF(x)
        sgn = 1 if err > 0 else -1
        errs.append((sgn, err, x))
    
    # Step1 - Alternating conditon.
    # 동일한 부호가 발생하는 경우 최대오차만 남기고 제거.
    res_step1 = [errs[0]]
    for i, err_data in enumerate(errs):
        if i == 0:
            continue
        
        compare_data = res_step1[-1]
        if compare_data[0] != err_data[0]:
            res_step1.append(err_data)
        else:
            if abs(compare_data[1]) < abs(err_data[1]):
                res_step1[-1] = err_data
        
    if len(res_step1) < sample_number:
        return [], []  
    
    # Step2 - Maximizing summation of errors
    res_step2 = [(d[0], abs(d[1]), d[2]) for d in res_step1]
    err_sum_neigh = []
    for i in range(0, len(res_step2)-1):
        err_sum_neigh.append((i, res_step2[i][1] + res_step2[i+1][1]))
    err_sum_neigh.sort(key=lambda t: t[1])
    
    while len(res_step2) > sample_number:
        if len(res_step2) == sample_number + 1:
            rm_index = 0 if res_step2[0][1] < res_step2[-1][1] else -1
            del res_step2[rm_index]
        
        elif len(res_step2) >= sample_number + 2:
            err_sum_neigh.append((-1, res_step2[0][1] + res_step2[-1][1]))
            err_sum_neigh.sort(key=lambda t: t[1])
            rm_index = err_sum_neigh[0][0]
            del res_step2[rm_index+1]
            del res_step2[rm_index]
    
    mpx = [data[2] for data in res_step2]
    mpy = [data[0] * data[1] for data in res_step2]
    return mpx, mpy

# 종료조건 판단
def decide_exit(errors: list, threshold: float, print_mode: str) -> bool:
    max_err, min_err = max(errors), min(errors)
    debug_print(f"Min/Max error\t\t: {min_err:.5f}, {max_err:.5f}", print_mode)
    if min_err != 0:
        exit_err = (max_err - min_err) / min_err
    else:
        exit_err = 0.0
    debug_print(f"Exit error\t\t: {exit_err}", print_mode)

    if exit_err <= threshold:
        return True
    else:
        return False
