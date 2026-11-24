# Write
import os, time

# Save
import pickle, shutil, tempfile
from pathlib import Path
from common.paths import OPS_TEMP

REMEZ_TEMP_PREFIX = "remez_step2_"

def cleanup_stale_remez_tempdirs() -> int:
    """
    이전 실행에서 비정상 종료로 남은 Remez 임시파일 및
    임시 디렉터리를 모두 제거한다.

    Returns:
        int: 제거한 임시 항목 개수
    """
    temp_root = OPS_TEMP
    temp_root.mkdir(parents=True, exist_ok=True)
    removed_count = 0

    try:
        candidates = list(temp_root.iterdir())
    except OSError as e:
        print(f"[Temp cleanup] Failed to scan '{temp_root}': {e}")
        return 0

    for path in candidates:
        if not path.name.startswith(REMEZ_TEMP_PREFIX):
            continue

        try:
            # 심볼릭 링크는 대상 경로를 따라가지 않고 링크만 제거한다.
            if path.is_symlink() or path.is_file():
                path.unlink()
            elif path.is_dir():
                shutil.rmtree(path)
            else:
                continue

            removed_count += 1

        except FileNotFoundError:
            # 다른 정리 작업에서 이미 삭제된 경우
            continue
        except OSError as e:
            print(f"[Temp cleanup] Failed to remove '{path}': {e}")

    return removed_count


def create_remez_tempdir() -> str:
    """
    Remez 디스크 저장용 임시 디렉터리를 생성한다.
    """
    OPS_TEMP.mkdir(parents=True, exist_ok=True)
    return tempfile.mkdtemp(prefix=REMEZ_TEMP_PREFIX, dir=OPS_TEMP)

class TqdmOverwriteFile:
    """
    파일에 진행바 1줄을 '덮어쓰기'로 유지하면서,
    일반 print 출력은 항상 파일 끝에 append 되도록 보장.
    """
    def __init__(self, fp, prefix="[PROGRESS] "):
        self.fp = fp
        self.prefix = prefix
        self._pos = None
        self._last_len = 0

    def write(self, s: str):
        if not s:
            return

        s = s.replace("\r", "").rstrip("\n")
        if not s:
            return

        # 최초: 추가 개행 없이 현재 위치를 진행바 라인 시작점으로 고정
        if self._pos is None:
            self._pos = self.fp.tell()

        # 진행바 라인 위치로 이동해서 덮어쓰기
        self.fp.seek(self._pos)
        line = self.prefix + s
        pad = max(0, self._last_len - len(line))
        self.fp.write(line + (" " * pad) + "\n")
        self.fp.flush()

        self._last_len = len(line)

        # 매우 중요: 이후 일반 print가 덮어쓰지 않도록 항상 EOF로 복귀
        self.fp.seek(0, os.SEEK_END)

    def flush(self):
        self.fp.flush()

def print_step_result(desc, current_n, total_n, start_time, calculated_count, completed_count, removed_len):
    step_dt = time.perf_counter() - start_time
    print(f"{desc}: DONE ({current_n}/{total_n}) | step_time={step_dt:.2f}s")
    print(f"\t{calculated_count + completed_count} routes newly calculated.")
    print(f"\t{completed_count} routes additionally completed.")
    # print(f"\t{removed_len} routes pruned.")
    print('-'*30)
    


class DiskRouteSet:
    def __init__(self, dir_path: str, count: int):
        self.dir_path = dir_path
        self.count = count

    def __len__(self):
        return self.count

    def chunk_files(self):
        base = Path(self.dir_path)
        return sorted(base.glob("chunk_*.pkl"))

    def cleanup(self):
        if os.path.exists(self.dir_path):
            shutil.rmtree(self.dir_path, ignore_errors=True)


def _save_routes_to_disk(routes: list, dir_path: str, chunk_size: int) -> DiskRouteSet:
    os.makedirs(dir_path, exist_ok=True)

    total = len(routes)
    file_idx = 0
    for i in range(0, total, chunk_size):
        chunk = routes[i:i + chunk_size]
        chunk_path = os.path.join(dir_path, f"chunk_{file_idx:06d}.pkl")
        with open(chunk_path, "wb") as f:
            pickle.dump(chunk, f, protocol=pickle.HIGHEST_PROTOCOL)
        file_idx += 1

    # 메모리 해제 유도
    routes.clear()
    return DiskRouteSet(dir_path, total)


def _iter_disk_routes(route_set: DiskRouteSet):
    for chunk_path in route_set.chunk_files():
        with open(chunk_path, "rb") as f:
            chunk = pickle.load(f)
        yield chunk_path, chunk


class RouteChunkWriter:
    """
    결과 route를 메모리에 무한정 쌓지 않고, 일정량마다 디스크로 flush.
    """
    def __init__(self, dir_path: str, chunk_size: int):
        self.dir_path = dir_path
        self.chunk_size = chunk_size
        self.buffer = []
        self.count = 0
        self.file_idx = 0
        os.makedirs(dir_path, exist_ok=True)

    def add(self, route):
        self.buffer.append(route)
        self.count += 1
        if len(self.buffer) >= self.chunk_size:
            self.flush()

    def extend(self, routes):
        for route in routes:
            self.add(route)

    def flush(self):
        if not self.buffer:
            return
        chunk_path = os.path.join(self.dir_path, f"chunk_{self.file_idx:06d}.pkl")
        with open(chunk_path, "wb") as f:
            pickle.dump(self.buffer, f, protocol=pickle.HIGHEST_PROTOCOL)
        self.buffer.clear()
        self.file_idx += 1

    def close(self) -> DiskRouteSet:
        self.flush()
        return DiskRouteSet(self.dir_path, self.count)
