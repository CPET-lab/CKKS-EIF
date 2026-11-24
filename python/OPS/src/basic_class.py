# basic_class.py
from OPS.src.basic import cal_complexity, detect_coeff_type
from OPD.src.basic_class import Complexity
from math import log2, ceil

class CalData:
    def __init__(self, title="", iter=0):
        self.title = title
        self.iter = iter
        self.depth = 0
        self.add_count = 0
        self.pmult = 0
        self.cmult = 0

    def add(self, c2: 'CalData', mode="add"):
        if mode == "add":
            self.depth += c2.depth
        elif mode == "compare":
            self.depth = self.depth if self.depth >= c2.depth else c2.depth
        self.add_count += c2.add_count
        self.pmult += c2.pmult
        self.cmult += c2.cmult
    
    def update(self, cmp: Complexity):
        self.depth = cmp.depth
        self.cmult = cmp.cmult
        self.pmult = cmp.pmult
        self.add_count = cmp.add
        
    def key(self, criteria="depth"):
        orders = {
            "depth": ("depth", "cmult", "pmult", "add_count"),
            "cmult": ("cmult", "depth", "pmult", "add_count"),
            "prune": ("cmult", "pmult", "add_count")
        }

        if criteria not in orders:
            raise ValueError(f"Unknown criteria: {criteria}")

        return tuple(getattr(self, name) for name in orders[criteria])

    def __lt__(self, other):
        return self.key("depth") < other.key("depth")

    def __eq__(self, other):
        return self.key("depth") == other.key("depth")

    def __le__(self, other):
        return self < other or self == other
    # def __lt__(self, other):
    #     return (self.depth, self.cmult, self.pmult, self.add_count) < (other.depth, other.cmult, other.pmult, other.add_count)
    
    # def __eq__(self, other):
    #     return (self.depth, self.cmult, self.pmult, self.add_count) == (other.depth, other.cmult, other.pmult, other.add_count)
    # def __le__(self, other):
    #     return self < other or self == other
    
    def print_params(self, title=False, iter=False):
        if title:
            print(f"Title:\t\t{self.title}")
        if iter:
            print(f"Iter:\t\t{self.iter}")
        print(f"DCPA: {self.depth}|{self.cmult}|{self.pmult}|{self.add_count}")

    
class remezData:
    def __init__(self, interval):
        self.coeff_log: list = []
        self.degree_log: list[int] = []
        self.interval_log: list = [interval]
        self.CalData_log: list[CalData] = []
        self.pre_log: list[int | float] = []
        self.pre_sep_log: list[tuple[float, float]] = []
        self.total_CalData = CalData()
        self.iter_cleanse = 0
        
        self.coeff_type_log: list = []
        self.time_log: list[float] = []
        self.total_time = 0.0

    def update(self, coeff, pre, pre_sep: tuple, es, interval, coeff_type):
        self.coeff_log.append(coeff)
        self.degree_log.append(len(coeff)-1)
        cd = CalData()
        cd.update(cal_complexity(coeff, es))
        self.CalData_log.append(cd)
        self.total_CalData.add(cd, "add")
        self.pre_log.append(pre)
        self.pre_sep_log.append(pre_sep)
        self.interval_log.append(interval)
        self.coeff_type_log.append(coeff_type)

    def update_time(self, time_cache):
        remain_depth = self.total_CalData.depth
        self.total_time = 0.0
        for coeff, coeff_type in zip(self.coeff_log, self.coeff_type_log):
            deg = len(coeff)-1
            time = time_cache[remain_depth][deg][coeff_type]
            self.time_log.append(time)
            self.total_time += time
            remain_depth -= ceil(log2(deg+1))
        
        
    def __lt__(self, other) -> bool:
        return (-self.pre_log[-1], self.total_CalData, self.pre_sep_log[-1][0], self.pre_sep_log[-1][1]) \
             < (-other.pre_log[-1], other.total_CalData, other.pre_sep_log[-1][0], other.pre_sep_log[-1][1])

    def __eq__(self, other) -> bool:
        return (self.total_CalData, self.pre_log[-1], self.pre_sep_log[-1][0], self.pre_sep_log[-1][1]) \
            == (other.total_CalData, other.pre_log[-1], other.pre_sep_log[-1][0], other.pre_sep_log[-1][1])

    def print_params(self):
        # print(f"precision: ", end='')
        # print(' -> '.join(map(str, self.pre_log)))
        self.total_CalData.print_params()
        print(f"coeff: ", end='')
        print(' -> '.join(map(str, [len(s)-1 for s in self.coeff_log])))
        print(f"Estimated time: {self.total_time}s.")
        # print(f"total complexity:")


def remez_key_pre(rd: remezData, criteria="depth"):
    return (
        -rd.pre_log[-1],
        rd.total_CalData.key(criteria),
        rd.pre_sep_log[-1][0],
        rd.pre_sep_log[-1][1],
    )

def remez_key_comp(rd: remezData, criteria="depth"):
    return (
        rd.total_CalData.key(criteria),
        -rd.pre_log[-1],
        rd.pre_sep_log[-1][0],
        rd.pre_sep_log[-1][1],
    )

def remez_key_coeff(rd: remezData, criteria="depth"):
    return (
        -rd.pre_log[-1],
        rd.total_CalData.key(criteria),
        tuple(len(coeff) for coeff in rd.coeff_log),
    )

# depth가 동일한 두 데이터를 비교할 때 사용
def remez_key_new(rd: remezData):
    return (
        -rd.pre_log[-1],
        tuple(cal_data.key("depth") for cal_data in rd.CalData_log),
    )
    
# pruning용 - depth가 동일한 두 데이터를 비교할 때
def remez_key_prune(rd: remezData):
    return (
        tuple(cal_data.key("prune") for cal_data in rd.CalData_log),
        -rd.pre_log[-1],
    )
