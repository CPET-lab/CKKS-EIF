import json
from pathlib import Path

def parse_time_cache(file_path: str | Path) -> dict:
    path = Path(file_path)

    if not path.is_file():
        raise FileNotFoundError(f"파일을 찾을 수 없습니다: {path}")

    try:
        with path.open("r", encoding="utf-8") as file:
            raw_data = json.load(file)
    except json.JSONDecodeError as error:
        raise ValueError(
            f"유효하지 않은 JSON 파일입니다: "
            f"line={error.lineno}, column={error.colno}"
        ) from error

    if not isinstance(raw_data, dict):
        raise TypeError("JSON의 최상위 자료형은 dict여야 합니다.")

    parsed_data: dict = {}

    for raw_first_key, raw_degree_dict in raw_data.items():
        try:
            first_key = int(raw_first_key)
        except (TypeError, ValueError) as error:
            raise ValueError(
                f"첫 번째 키를 정수로 변환할 수 없습니다: {raw_first_key!r}"
            ) from error

        if not isinstance(raw_degree_dict, dict):
            raise TypeError(
                f"첫 번째 키 {first_key}의 값은 dict여야 합니다."
            )

        parsed_data[first_key] = {}

        for raw_degree, raw_type_dict in raw_degree_dict.items():
            try:
                degree = int(raw_degree)
            except (TypeError, ValueError) as error:
                raise ValueError(
                    f"두 번째 키를 정수로 변환할 수 없습니다: "
                    f"first_key={first_key}, degree={raw_degree!r}"
                ) from error

            if not isinstance(raw_type_dict, dict):
                raise TypeError(
                    f"({first_key}, {degree})의 값은 dict여야 합니다."
                )

            parsed_data[first_key][degree] = {}

            for polynomial_type, raw_time in raw_type_dict.items():
                if polynomial_type not in {"all", "even", "odd"}:
                    raise ValueError(
                        f"지원하지 않는 다항식 유형입니다: "
                        f"first_key={first_key}, degree={degree}, "
                        f"type={polynomial_type!r}"
                    )

                # bool은 int의 하위 클래스이므로 별도로 제외한다.
                if isinstance(raw_time, bool) or not isinstance(
                    raw_time, (int, float)
                ):
                    raise TypeError(
                        f"실행시간은 숫자여야 합니다: "
                        f"first_key={first_key}, degree={degree}, "
                        f"type={polynomial_type}, value={raw_time!r}"
                    )

                parsed_data[first_key][degree][polynomial_type] = float(raw_time)

    # JSON상의 문자열 정렬이 아니라 정수 기준으로 정렬한다.
    return {
        first_key: {
            degree: parsed_data[first_key][degree]
            for degree in sorted(parsed_data[first_key])
        }
        for first_key in sorted(parsed_data)
    }


def validate_time_cache(time_cache: dict) -> None:
    """
    각 차수에 필요한 유형이 올바르게 들어 있는지 검사한다.

    짝수 차수: {"all", "even"}
    홀수 차수: {"all", "odd"}
    """
    for first_key, degree_dict in time_cache.items():
        for degree, type_dict in degree_dict.items():
            expected_types = {
                "all",
                "even" if degree % 2 == 0 else "odd",
            }
            actual_types = set(type_dict)

            if actual_types != expected_types:
                raise ValueError(
                    f"유형 구성이 올바르지 않습니다: "
                    f"first_key={first_key}, degree={degree}, "
                    f"expected={sorted(expected_types)}, "
                    f"actual={sorted(actual_types)}"
                )


