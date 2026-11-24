"""OPD 탐색 캐시와 JSON 저장 형식의 자료형 명세.

이 모듈은 두 종류의 자료형을 구분한다.

1. 메모리 자료형
   OPD 탐색 중 사용하는 SearchCacheKey, SearchResult와 power 상태다.
   NamedTuple을 사용하므로 기존 일반 tuple과 동일하게 해시, 비교, unpack이
   가능하다.
2. 직렬화 자료형
   Serialized* TypedDict는 :mod:`save`가 JSON으로 기록하는 dictionary의
   필수 key와 value 형태를 명시한다. TypedDict는 실행 시 dictionary 그대로다.
"""

from typing import (
    TYPE_CHECKING,
    Any,
    Literal,
    NamedTuple,
    NotRequired,
    TypeAlias,
    TypedDict,
)

if TYPE_CHECKING:
    from .basic_class import Complexity, Decomp
else:
    # get_type_hints()도 전방 참조를 해석할 수 있게 런타임 이름을 제공한다.
    # 실제 캐시 값의 타입 제한은 정적 분석 시 위 import가 담당한다.
    Complexity = Any
    Decomp = Any


# ---------------------------------------------------------------------------
# OPD 실행 중 여러 모듈이 공유하는 원자 자료형
# ---------------------------------------------------------------------------

# 다항식 계수는 상수항부터 최고차항 순서로 저장된다.
CoefficientType: TypeAlias = Literal["0", "I", "F"]
NumericCoefficient: TypeAlias = int | float
PolynomialType: TypeAlias = list[CoefficientType]
NumericPolynomial: TypeAlias = list[NumericCoefficient]

# power는 지수, depth는 해당 순수 power x^power가 실제로 생성된 depth다.
Power: TypeAlias = int
PowerDepth: TypeAlias = int
PowerDepthPair: TypeAlias = tuple[Power, PowerDepth]
PowerDepths: TypeAlias = dict[Power, PowerDepth]

# 캐시 key에 넣을 때 dictionary를 hashable한 frozenset으로 바꾼 표현이다.
PowerState: TypeAlias = frozenset[PowerDepthPair]

# pure: x^lhs * x^rhs, coeff: ax^lhs * x^rhs,
# coeff_direct: a * x^rhs를 의미한다.
RouteKind: TypeAlias = Literal["pure", "coeff", "coeff_direct"]
Route: TypeAlias = tuple[int, int]
RouteOp: TypeAlias = tuple[RouteKind, int, int]
EvalOrder: TypeAlias = Literal[
    "terminal",
    "p_only",
    "q_only",
    "p_then_q",
    "q_then_p",
]

# 분해하지 않는 terminal 다항식의 비영 항 하나를 계산하는 메모리 평가계획
class TermPlan(TypedDict):
    degree: int
    coeff_type: CoefficientType
    multA: bool
    route_ops: list[RouteOp]


# ---------------------------------------------------------------------------
# 메모리 탐색 캐시
# ---------------------------------------------------------------------------

# 다항식 키 -> 타입 또는 값 tuple.
PolynomialKey: TypeAlias = tuple[CoefficientType, ...] | tuple[float, ...]

# 같은 OPD 부분문제를 식별하는 hashable key.
class SearchCacheKey(NamedTuple):
    """key, depth 및 option이 동일할 경우 같은 것으로 판단함.
        
    poly_key는 타입 다항식이면 0/I/F tuple, 수치 다항식이면 실제 계수
    tuple이다. 입력 power 상태는 별도의 :class:`SearchResult`에 저장한다.
    """

    poly_key: PolynomialKey
    is_type_poly: bool
    divAssumption: bool
    depth_budget: int

# 특정 입력 power 상태에서 얻은 decomposition 후보 클래스.
class SearchResult(NamedTuple):
    """
    decomposition과 complexity는 순환 import를 피하기 위해 전방 참조로
    표기한다. 실제 값은 각각 basic_class.Decomp와 Complexity다.
    """

    power_state: PowerState
    decomposition: "Decomp"
    complexity: "Complexity"


# 캐시 - (키, list(후보))의 형태.
SearchCache: TypeAlias = dict[SearchCacheKey, list[SearchResult]]


# ---------------------------------------------------------------------------
# JSON 직렬화 dictionary
# ---------------------------------------------------------------------------

# RouteOp의 JSON 객체 표현
class SerializedRouteOp(TypedDict):
    kind: RouteKind
    lhs: int
    rhs: int

# TermPlan의 JSON 객체 표현
class SerializedTermPlan(TypedDict):
    degree: int
    coeff_type: CoefficientType
    multA: bool
    route_ops: list[SerializedRouteOp]

# Complexity의 JSON 객체 표현
class SerializedComplexity(TypedDict):
    depth: int
    cmult: int
    pmult: int
    add: int

# XI의 JSON 객체 표현
class SerializedXI(TypedDict):
    n: int
    multA: bool
    route: list[list[int]]
    route_ops: list[SerializedRouteOp]
    power_depths: list[list[int]]
    depth: int
    cmult: int
    pmult: int
    # 이전 저장 포맷을 읽을 때만 사용되는 선택 필드다.
    made_powers: NotRequired[list[int]]

# 재귀 decomposition 트리 노드의 JSON 객체 표현
class SerializedDecomp(TypedDict):
    is_leaf: bool
    coeff_type: list[CoefficientType]
    eval_order: EvalOrder
    complexity: SerializedComplexity
    xi: SerializedXI | None
    term_plans: list[SerializedTermPlan]
    p: "SerializedDecomp | None"
    q: "SerializedDecomp | None"

# SearchCacheKey의 JSON 객체 표현
class SerializedSearchKey(TypedDict):
    poly_key: list[str | float]
    is_type_poly: bool
    divAssumption: bool
    depth_budget: int

# SearchResult의 JSON 객체 표현
class SerializedSearchResult(TypedDict):
    input_power_depths: list[list[int]]
    complexity: SerializedComplexity
    decomposition: SerializedDecomp
    # 이전 저장 포맷을 읽을 때만 사용되는 선택 필드다.
    input_powers: NotRequired[list[int]]

# 하나의 SearchCacheKey와 그 후보 목록을 묶은 JSON entry
class SerializedSearchEntry(TypedDict):
    key: SerializedSearchKey
    results: list[SerializedSearchResult]

# save_opd_caches가 파일 최상위에 기록하는 구조
class SerializedOPDCaches(TypedDict):
    initial: list[SerializedSearchEntry]
    intermediate: list[SerializedSearchEntry]


# 계수타입 문자열(예: F0IF) 하나에 최적 분해 하나를 저장하는 단순 캐시다.
DecompositionCache: TypeAlias = dict[str, SerializedDecomp]


__all__ = [
    "CoefficientType",
    "DecompositionCache",
    "EvalOrder",
    "NumericCoefficient",
    "NumericPolynomial",
    "PolynomialKey",
    "PolynomialType",
    "PowerDepthPair",
    "PowerDepths",
    "PowerState",
    "Route",
    "RouteKind",
    "RouteOp",
    "SearchCache",
    "SearchCacheKey",
    "SearchResult",
    "SerializedComplexity",
    "SerializedDecomp",
    "SerializedOPDCaches",
    "SerializedRouteOp",
    "SerializedSearchEntry",
    "SerializedSearchKey",
    "SerializedSearchResult",
    "SerializedTermPlan",
    "SerializedXI",
    "TermPlan",
]
