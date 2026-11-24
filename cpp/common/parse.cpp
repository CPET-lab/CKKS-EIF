// parse.cpp
#include "parse.h"

static double stod_safe(const std::string& s) {
    // 필요시 안전 변환용 헬퍼
    return std::stod(s);
}

std::map<std::string, fs::path> findFilesByPrefix(
    const fs::path& directory,
    std::string_view prefix)
{
    const fs::path targetDirectory =
        directory.is_absolute()
            ? directory
            : fs::current_path() / directory;

    std::map<std::string, fs::path> result;

    if (!fs::exists(targetDirectory) ||
        !fs::is_directory(targetDirectory)) {
        return result;
    }

    for (const auto& entry :
         fs::directory_iterator(targetDirectory)) {
        if (!entry.is_regular_file()) {
            continue;
        }

        const auto text =
            extractLastString(entry.path(), prefix);

        if (!text) {
            continue;
        }

        result.emplace(*text, entry.path());
    }

    return result;
}

std::optional<std::string> extractLastString(
    const fs::path& path,
    std::string_view prefix)
{
    const std::string stem = path.stem().string();

    if (!stem.starts_with(prefix)) {
        return std::nullopt;
    }

    const std::string text = stem.substr(prefix.size());

    if (text.empty()) {
        return std::nullopt;
    }

    return text;
}

std::string remez_coeff_filename(int p_num, int e_num, const std::string& criteria)
{
    if(criteria != "depth" && criteria != "time")
        throw std::invalid_argument("criteria must be depth or time");
    return "coeff_" + std::to_string(p_num) + "_" + std::to_string(e_num)
        + "_" + criteria + ".txt";
}

std::vector<std::vector<double>> parse_remez_coeff(
    int p_num, int e_num, const std::string& criteria, bool print, const std::string& dir)
{
    return parse_remez_coeff(p_num, e_num, print, dir,
        remez_coeff_filename(p_num, e_num, criteria));
}

std::vector<std::vector<double>> parse_remez_coeff(
    int p_num, int e_num, bool print, std::string dir, std::string filename)
{
    if(filename != remez_coeff_filename(p_num, e_num, "depth")
        && filename != remez_coeff_filename(p_num, e_num, "time"))
        throw std::invalid_argument("Unexpected coefficient filename: " + filename);

    int line_count = 0;  // 전체 줄 개수

    fs::path file_path = fs::path(dir) / filename;

    // 2) 계수 파일 파싱
    std::ifstream infile(file_path);
    if (!infile)
        throw std::runtime_error("파일을 열 수 없습니다: " + fs::absolute(file_path).string());

    std::vector<std::vector<double>> result;
    std::string line;

    // 계수 추출 패턴: (계수)x^(지수)
    const std::regex term_re(
        R"(([+\-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+\-]?\d+)?))" // coeff (그룹1)
        R"(x\^)"                                              // x^
        R"((\d+))"                                            // exponent (그룹2)
    );

    // s 계산용 패턴: "+3x^2-2x^3" (공백 가변 허용, 앞 '+'는 선택)
    const std::regex s_re(R"((?:\+)?\s*3\s*x\^2\s*-\s*2\s*x\^3)");

    // 파일 순회
    while (std::getline(infile, line)) {
        ++line_count;  // 줄 개수 증가

        // (exp -> coeff) 맵 구성
        std::map<int, double> exp_to_coeff;
        int max_exp = -1;

        for (std::sregex_iterator it(line.begin(), line.end(), term_re);
            it != std::sregex_iterator(); ++it)
        {
            const std::string coeff_str = (*it)[1].str();
            const std::string exp_str   = (*it)[2].str();

            double c = stod_safe(coeff_str);
            int    e = std::stoi(exp_str);

            exp_to_coeff[e] = c;
            if (e > max_exp) max_exp = e;
        }

        // 항이 하나도 없으면 빈 벡터
        if (max_exp < 0) {
            result.emplace_back();
            continue;
        }

        // 누락 차수 0.0으로 채움
        std::vector<double> coeffs(static_cast<size_t>(max_exp + 1), 0.0);
        for (const auto& kv : exp_to_coeff) {
            coeffs[static_cast<size_t>(kv.first)] = kv.second;
        }

        result.push_back(std::move(coeffs));
    }
    return result;
}

int calcuate_remez_depth(std::vector<std::vector<double>> coeffs)
{
    int total = 0;
    for(auto coeff: coeffs)
    {
        total += ceil(log2(coeff.size()));
    }
    return total;
}

nlohmann::json load_decomp_cache(const std::string& filename)
{
    std::ifstream f(filename);
    
    // 파일 열기 실패 확인
    if (!f.is_open()) {
        std::cerr << "Failed to open " << filename << std::endl;
        return nullptr; // 실패 시 null 상태의 json 반환
    }

    nlohmann::json cache_root;
    try {
        cache_root = nlohmann::json::parse(f);
    } catch (nlohmann::json::parse_error& e) {
        std::cerr << "JSON Parse Error: " << e.what() << std::endl;
        return nullptr;
    }

    return cache_root;
}

std::string get_poly_type_key(const std::vector<double>& coeff)
{
    std::string key;
    key.reserve(coeff.size());

    for (double v : coeff) {
        double ip;
        if(v == 0.0)
            key.push_back('0');
        else if(std::modf(v, &ip) == 0.0)
            key.push_back('I');
        else
            key.push_back('F');

        // key.push_back(v == 0.0 ? '0' : 'F');
    }
    return key;
}

/*
 * Legacy decomposition-cache parser.
 *
 * 이 구현은 is_leaf/xi/p/q와 의미가 구분되지 않은 xi.route만 저장하던
 * 과거 cache 전용이다. 현재 cache의 leaf 실행계획, typed route 및 평가
 * 순서를 보존할 수 없으므로 실행 코드에서는 사용하지 않는다.
 *
std::shared_ptr<Decomp> reconstruct_decomp_from_cache(
    const json& plan,
    Poly poly,
    int depth)
{
    if (plan["is_leaf"].get<bool>())
        return std::make_shared<Decomp>(poly.coeff, Complexity());

    json plan_xi = plan["xi"];
    XI xi = XI(plan_xi.value("multA", false), plan_xi.value("n", 0));

    std::set<int> made_powers = {0, 1};
    std::vector<std::pair<int, int>> routes;
    for (auto& route : plan_xi.value("route", json::array()))
        routes.push_back({route[0], route[1]});

    for (std::size_t i = 0; i < routes.size(); ++i) {
        if (i + 1 == routes.size() && xi.multA)
            break;
        made_powers.insert(routes[i].first + routes[i].second);
    }

    xi.coeff = poly.coeff.back();
    xi.add_routes(routes, made_powers);

    std::pair<Poly, Poly> polys = poly.seperate_poly(xi.n, xi.multA);
    std::shared_ptr<Decomp> dcmp_p = reconstruct_decomp_from_cache(
        plan["p"],
        polys.first,
        depth + 1
    );
    std::shared_ptr<Decomp> dcmp_q;
    if (!polys.second.coeff.empty()) {
        dcmp_q = reconstruct_decomp_from_cache(
            plan["q"],
            polys.second,
            depth + 1
        );
    }

    std::shared_ptr<Decomp> result = std::make_shared<Decomp>(
        poly.coeff,
        Complexity()
    );
    result->total_routes = xi.route;
    result->update(xi, dcmp_p, dcmp_q);
    result->made_powers.insert(
        dcmp_p->made_powers.begin(),
        dcmp_p->made_powers.end()
    );
    result->merge_route(dcmp_p->total_routes);
    if (dcmp_q) {
        result->made_powers.insert(
            dcmp_q->made_powers.begin(),
            dcmp_q->made_powers.end()
        );
        result->merge_route(dcmp_q->total_routes);
    }
    return result;
}
*/

namespace
{

struct ParsedXIData
{
    bool exists = false;
    bool multA = false;
    int n = 0;
    int depth = 0;
    int cmult = 0;
    int pmult = 0;
    std::vector<std::pair<int, int>> routes;
    std::vector<RouteOp> route_ops;
    std::set<int> made_powers = {0, 1};
};

std::runtime_error make_plan_error(int depth, const std::string& message)
{
    return std::runtime_error(
        "Invalid decomposition plan at depth "
        + std::to_string(depth)
        + ": "
        + message
    );
}

RouteKind parse_route_kind(
    const std::string& kind,
    int depth,
    const std::string& location
)
{
    if (kind == "pure")
        return RouteKind::Pure;
    if (kind == "coeff")
        return RouteKind::Coeff;
    if (kind == "coeff_direct")
        return RouteKind::CoeffDirect;

    throw make_plan_error(
        depth,
        location + " has invalid route kind: " + kind
    );
}

std::vector<RouteOp> parse_route_ops(
    const json& owner,
    bool legacy_multA,
    int depth,
    const std::string& location
)
{
    std::vector<RouteOp> result;

    if (owner.contains("route_ops") && !owner.at("route_ops").is_null()) {
        const json& raw_ops = owner.at("route_ops");
        if (!raw_ops.is_array())
            throw make_plan_error(depth, location + ".route_ops is not an array");

        result.reserve(raw_ops.size());
        for (std::size_t i = 0; i < raw_ops.size(); ++i) {
            const json& raw = raw_ops.at(i);
            if (!raw.is_object()) {
                throw make_plan_error(
                    depth,
                    location + ".route_ops[" + std::to_string(i) + "] is not an object"
                );
            }

            if (!raw.contains("lhs") || !raw.at("lhs").is_number_integer()
                || !raw.contains("rhs") || !raw.at("rhs").is_number_integer()) {
                throw make_plan_error(
                    depth,
                    location + ".route_ops[" + std::to_string(i)
                    + "] must contain integer lhs/rhs"
                );
            }

            const int lhs = raw.at("lhs").get<int>();
            const int rhs = raw.at("rhs").get<int>();
            if (lhs < 0 || rhs < 0) {
                throw make_plan_error(
                    depth,
                    location + ".route_ops[" + std::to_string(i)
                    + "] contains a negative power"
                );
            }

            result.push_back(RouteOp{
                parse_route_kind(
                    raw.value("kind", std::string{}),
                    depth,
                    location + ".route_ops[" + std::to_string(i) + "]"
                ),
                lhs,
                rhs
            });
        }
        return result;
    }

    // schema v1 compatibility
    if (!owner.contains("route") || owner.at("route").is_null())
        return result;

    const json& raw_routes = owner.at("route");
    if (!raw_routes.is_array())
        throw make_plan_error(depth, location + ".route is not an array");

    result.reserve(raw_routes.size());
    for (std::size_t i = 0; i < raw_routes.size(); ++i) {
        const json& raw = raw_routes.at(i);
        if (!raw.is_array() || raw.size() != 2
            || !raw.at(0).is_number_integer()
            || !raw.at(1).is_number_integer()) {
            throw make_plan_error(
                depth,
                location + ".route[" + std::to_string(i)
                + "] must contain two integer powers"
            );
        }

        RouteKind kind = RouteKind::Pure;
        if (legacy_multA && i + 1 == raw_routes.size())
            kind = RouteKind::Coeff;

        result.push_back(RouteOp{
            kind,
            raw.at(0).get<int>(),
            raw.at(1).get<int>()
        });
    }
    return result;
}

void apply_route_ops(
    const std::vector<RouteOp>& route_ops,
    std::set<int>& available_powers,
    int depth,
    const std::string& location
)
{
    std::set<int> coefficient_powers;

    for (std::size_t i = 0; i < route_ops.size(); ++i) {
        const RouteOp& op = route_ops[i];
        const std::string op_location = location + "[" + std::to_string(i) + "]";

        if (op.kind == RouteKind::Pure) {
            if (!available_powers.contains(op.lhs)) {
                throw make_plan_error(
                    depth,
                    op_location + " uses unavailable pure power x^"
                    + std::to_string(op.lhs)
                );
            }
            if (!available_powers.contains(op.rhs)) {
                throw make_plan_error(
                    depth,
                    op_location + " uses unavailable pure power x^"
                    + std::to_string(op.rhs)
                );
            }
            available_powers.insert(op.lhs + op.rhs);
            continue;
        }

        if (op.kind == RouteKind::CoeffDirect) {
            if (!available_powers.contains(op.rhs)) {
                throw make_plan_error(
                    depth,
                    op_location + " uses unavailable pure power x^"
                    + std::to_string(op.rhs)
                );
            }
            coefficient_powers.insert(op.rhs);
            continue;
        }

        // (a*x^lhs) * x^rhs
        if (!available_powers.contains(op.rhs)) {
            throw make_plan_error(
                depth,
                op_location + " uses unavailable pure power x^"
                + std::to_string(op.rhs)
            );
        }
        if (!coefficient_powers.contains(op.lhs)
            && !available_powers.contains(op.lhs)) {
            throw make_plan_error(
                depth,
                op_location + " uses unavailable coefficient base x^"
                + std::to_string(op.lhs)
            );
        }
        coefficient_powers.insert(op.lhs + op.rhs);
    }
}

Complexity parse_complexity(const json& plan, int depth)
{
    Complexity comp;
    if (!plan.contains("complexity") || plan.at("complexity").is_null())
        return comp;

    const json& raw = plan.at("complexity");
    if (!raw.is_object())
        throw make_plan_error(depth, "complexity is not an object");

    comp.insert_value(
        raw.value("depth", 0),
        raw.value("cmult", 0),
        raw.value("pmult", 0),
        raw.value("add", 0)
    );
    return comp;
}

ParsedXIData parse_xi_data(
    const json& plan,
    const std::set<int>& incoming_powers,
    int depth,
    bool required
)
{
    ParsedXIData result;
    result.made_powers = incoming_powers;

    if (!plan.contains("xi") || plan.at("xi").is_null()) {
        if (required)
            throw make_plan_error(depth, "missing xi");
        return result;
    }

    const json& raw = plan.at("xi");
    if (!raw.is_object())
        throw make_plan_error(depth, "xi is not an object");

    result.exists = true;
    result.multA = raw.value("multA", false);
    result.n = raw.value("n", 0);
    result.depth = raw.value("depth", 0);
    result.cmult = raw.value("cmult", 0);
    result.pmult = raw.value("pmult", 0);

    if (result.n < 0)
        throw make_plan_error(depth, "xi.n is negative");

    result.route_ops = parse_route_ops(raw, result.multA, depth, "xi");
    apply_route_ops(result.route_ops, result.made_powers, depth, "xi.route_ops");

    result.routes.reserve(result.route_ops.size());
    for (const RouteOp& op : result.route_ops)
        result.routes.emplace_back(op.lhs, op.rhs);

    if (raw.contains("made_powers") && !raw.at("made_powers").is_null()) {
        const json& stored = raw.at("made_powers");
        if (!stored.is_array())
            throw make_plan_error(depth, "xi.made_powers is not an array");

        for (const json& value : stored) {
            if (!value.is_number_integer())
                throw make_plan_error(depth, "xi.made_powers contains a non-integer");
            const int power = value.get<int>();
            if (power < 0)
                throw make_plan_error(depth, "xi.made_powers contains a negative power");
            result.made_powers.insert(power);
        }
    }

    return result;
}

std::vector<TermPlan> parse_term_plans(const json& plan, int depth)
{
    std::vector<TermPlan> result;
    if (!plan.contains("term_plans") || plan.at("term_plans").is_null())
        return result;

    const json& raw_plans = plan.at("term_plans");
    if (!raw_plans.is_array())
        throw make_plan_error(depth, "term_plans is not an array");

    result.reserve(raw_plans.size());
    for (std::size_t i = 0; i < raw_plans.size(); ++i) {
        const json& raw = raw_plans.at(i);
        if (!raw.is_object()) {
            throw make_plan_error(
                depth,
                "term_plans[" + std::to_string(i) + "] is not an object"
            );
        }
        if (!raw.contains("degree") || !raw.at("degree").is_number_integer()) {
            throw make_plan_error(
                depth,
                "term_plans[" + std::to_string(i) + "].degree is invalid"
            );
        }

        const int degree = raw.at("degree").get<int>();
        if (degree < 0)
            throw make_plan_error(depth, "term plan contains a negative degree");

        result.push_back(TermPlan{
            degree,
            raw.value("coeff_type", std::string{"0"}),
            raw.value("multA", false),
            parse_route_ops(
                raw,
                raw.value("multA", false),
                depth,
                "term_plans[" + std::to_string(i) + "]"
            )
        });
    }
    return result;
}

bool read_leaf_flag(const json& plan, int depth)
{
    if (!plan.contains("is_leaf") || !plan.at("is_leaf").is_boolean())
        throw make_plan_error(depth, "missing or invalid is_leaf");
    return plan.at("is_leaf").get<bool>();
}

std::string read_eval_order(const json& plan, bool is_leaf, int depth)
{
    if (!plan.contains("eval_order") || !plan.at("eval_order").is_string())
        throw make_plan_error(depth, "missing or invalid eval_order");

    const std::string order = plan.at("eval_order").get<std::string>();

    static const std::set<std::string> valid = {
        "terminal", "p_only", "q_only", "p_then_q", "q_then_p"
    };
    if (!valid.contains(order))
        throw make_plan_error(depth, "invalid eval_order: " + order);
    if (is_leaf && order != "terminal")
        throw make_plan_error(depth, "leaf node must use eval_order=terminal");
    if (!is_leaf && order == "terminal")
        throw make_plan_error(depth, "non-leaf node cannot use eval_order=terminal");
    return order;
}

void append_unique_routes(
    std::vector<std::pair<int, int>>& destination,
    const std::vector<std::pair<int, int>>& source
)
{
    std::set<std::pair<int, int>> existing(
        destination.begin(),
        destination.end()
    );
    for (const auto& route : source) {
        if (existing.insert(route).second)
            destination.push_back(route);
    }
}

void append_pure_routes(
    std::vector<std::pair<int, int>>& destination,
    const std::vector<RouteOp>& ops
)
{
    std::set<std::pair<int, int>> existing(
        destination.begin(),
        destination.end()
    );
    for (const RouteOp& op : ops) {
        if (op.kind != RouteKind::Pure)
            continue;
        const std::pair<int, int> route{op.lhs, op.rhs};
        if (existing.insert(route).second)
            destination.push_back(route);
    }
}

std::shared_ptr<Decomp> reconstruct_decomp_from_cache_impl(
    const json& plan,
    Poly poly,
    const std::set<int>& incoming_powers,
    int depth
)
{
    if (!plan.is_object())
        throw make_plan_error(depth, "plan is null or not an object");

    // 버전 및 OPD option metadata 없이 현재 실행계획 필드를 직접 읽는다.

    const bool is_leaf = read_leaf_flag(plan, depth);
    const std::string eval_order = read_eval_order(plan, is_leaf, depth);
    ParsedXIData xi_data = parse_xi_data(
        plan,
        incoming_powers,
        depth,
        true
    );
    const std::vector<TermPlan> term_plans = parse_term_plans(plan, depth);
    const Complexity comp = parse_complexity(plan, depth);

    // terminal term plan의 순수 power 생성도 최종 available power에 반영한다.
    if (is_leaf && !term_plans.empty()) {
        std::set<int> terminal_available = incoming_powers;
        for (std::size_t i = 0; i < term_plans.size(); ++i) {
            apply_route_ops(
                term_plans[i].route_ops,
                terminal_available,
                depth,
                "term_plans[" + std::to_string(i) + "].route_ops"
            );
        }
        terminal_available.insert(
            xi_data.made_powers.begin(),
            xi_data.made_powers.end()
        );
        xi_data.made_powers = std::move(terminal_available);
    }

    XI xi(xi_data.multA, xi_data.n);
    if (!poly.coeff.empty())
        xi.coeff = poly.coeff.back();
    xi.add_routes(
        xi_data.routes,
        xi_data.made_powers,
        xi_data.route_ops
    );
    xi.depth = xi_data.depth;
    xi.cmult = xi_data.cmult;
    xi.pmult = xi_data.pmult;
    xi.add_count = xi.cmult;

    if (is_leaf) {
        auto result = std::make_shared<Decomp>(poly.coeff, comp, xi);
        result->has_execution_plan = true;
        result->eval_order = eval_order;
        result->term_plans = term_plans;
        result->made_powers = xi_data.made_powers;
        result->total_routes.clear();

        if (!term_plans.empty()) {
            for (const TermPlan& term : term_plans)
                append_pure_routes(result->total_routes, term.route_ops);
        } else {
            append_pure_routes(result->total_routes, xi_data.route_ops);
        }
        return result;
    }

    if (!xi_data.exists)
        throw make_plan_error(depth, "non-leaf node has no xi");
    if (xi_data.n == 0 && !xi_data.multA)
        throw make_plan_error(depth, "non-leaf node has n == 0 and multA == false");
    if (poly.coeff.empty())
        throw make_plan_error(depth, "non-leaf node received an empty polynomial");

    std::pair<Poly, Poly> split = poly.seperate_poly(xi.n, xi.multA);
    Poly poly_p = std::move(split.first);
    Poly poly_q = std::move(split.second);

    std::shared_ptr<Decomp> dcmp_p;
    std::shared_ptr<Decomp> dcmp_q;
    std::set<int> current_powers = xi_data.made_powers;

    const auto reconstruct_p = [&]() {
        if (poly_p.coeff.empty())
            return;
        if (!plan.contains("p") || plan.at("p").is_null())
            throw make_plan_error(depth, "p polynomial exists but p plan is null");

        dcmp_p = reconstruct_decomp_from_cache_impl(
            plan.at("p"),
            poly_p,
            current_powers,
            depth + 1
        );
        current_powers.insert(
            dcmp_p->made_powers.begin(),
            dcmp_p->made_powers.end()
        );
    };

    const auto reconstruct_q = [&]() {
        if (poly_q.coeff.empty())
            return;
        if (!plan.contains("q") || plan.at("q").is_null())
            throw make_plan_error(depth, "q polynomial exists but q plan is null");

        dcmp_q = reconstruct_decomp_from_cache_impl(
            plan.at("q"),
            poly_q,
            current_powers,
            depth + 1
        );
        current_powers.insert(
            dcmp_q->made_powers.begin(),
            dcmp_q->made_powers.end()
        );
    };

    if (eval_order == "q_then_p" || eval_order == "q_only") {
        reconstruct_q();
        if (eval_order != "q_only")
            reconstruct_p();
    } else {
        reconstruct_p();
        if (eval_order != "p_only")
            reconstruct_q();
    }

    auto result = std::make_shared<Decomp>(poly.coeff, comp, xi);
    result->update(xi, dcmp_p, dcmp_q);
    result->comp = comp;
    result->has_execution_plan = true;
    result->eval_order = eval_order;
    result->term_plans = term_plans;
    result->made_powers = current_powers;
    result->total_routes.clear();

    append_pure_routes(result->total_routes, xi_data.route_ops);
    if (eval_order == "q_then_p" || eval_order == "q_only") {
        if (dcmp_q)
            append_unique_routes(result->total_routes, dcmp_q->total_routes);
        if (dcmp_p)
            append_unique_routes(result->total_routes, dcmp_p->total_routes);
    } else {
        if (dcmp_p)
            append_unique_routes(result->total_routes, dcmp_p->total_routes);
        if (dcmp_q)
            append_unique_routes(result->total_routes, dcmp_q->total_routes);
    }

    return result;
}

} // namespace

std::shared_ptr<Decomp> reconstruct_decomp_from_cache(
    const json& plan,
    Poly poly,
    int depth
)
{
    return reconstruct_decomp_from_cache_impl(
        plan,
        std::move(poly),
        std::set<int>{0, 1},
        depth
    );
}

// Class EvalStep
EvalStep::EvalStep() = default;

EvalStep::EvalStep(std::shared_ptr<Decomp> dcmp)
{
    if (!dcmp)
        throw std::invalid_argument("EvalStep: dcmp is null");

    this->add_value('P', 1, "x");

    if (!dcmp->has_execution_plan) {
        for (const auto& route : dcmp->total_routes) {
            const int target = route.first + route.second;
            const std::string target_key = "P" + std::to_string(target);
            if (this->powers.contains(target_key))
                continue;

            this->add_step(
                "P" + std::to_string(route.first),
                "P" + std::to_string(route.second),
                'x'
            );
        }
    }

    serialize_dcmp(std::move(dcmp));
}

std::string EvalStep::serialize_dcmp(std::shared_ptr<Decomp> dcmp)
{
    if (!dcmp)
        throw std::invalid_argument("EvalStep::serialize_dcmp: dcmp is null");

    if (!dcmp->has_execution_plan) {
        XI xi = dcmp->xi;
        std::string coeff_key, xi_key, axi_key;

        if (!dcmp->dcmp_p && !dcmp->dcmp_q) {
            std::string res_key;
            for (int i = static_cast<int>(dcmp->coeff.size()) - 1; i >= 0; --i) {
                if (dcmp->coeff[i] == 0.0)
                    continue;

                std::string term_key;
                if (i == 0) {
                    term_key = this->add_value(
                        'C',
                        this->coeff_count++,
                        std::format("{}", dcmp->coeff[i])
                    );
                } else if (dcmp->coeff[i] == 1.0) {
                    term_key = "P" + std::to_string(i);
                } else {
                    const std::string ckey = this->add_value(
                        'C',
                        this->coeff_count++,
                        std::format("{}", dcmp->coeff[i])
                    );
                    term_key = this->add_step(
                        ckey,
                        "P" + std::to_string(i),
                        'x'
                    );
                }

                res_key = res_key.empty()
                    ? term_key
                    : this->add_step(res_key, term_key, '+');
            }
            return res_key;
        }

        if (xi.multA) {
            coeff_key = this->add_value(
                'C',
                this->coeff_count++,
                std::format("{}", xi.coeff)
            );
            xi_key = "P" + std::to_string(xi.n);
            axi_key = this->add_step(coeff_key, xi_key, 'x');
        } else {
            axi_key = "P" + std::to_string(xi.n);
        }

        std::string result_key = axi_key;
        if (dcmp->dcmp_p) {
            const std::string p_key = serialize_dcmp(dcmp->dcmp_p);
            result_key = this->add_step(axi_key, p_key, 'x');
        }
        if (dcmp->dcmp_q) {
            const std::string q_key = serialize_dcmp(dcmp->dcmp_q);
            result_key = this->add_step(result_key, q_key, '+');
        }
        return result_key;
    }

    const auto power_key = [](int exponent) {
        return "P" + std::to_string(exponent);
    };

    const auto has_power = [&](int exponent) {
        return this->powers.contains(power_key(exponent));
    };

    const auto require_power = [&](int exponent, const std::string& context) {
        if (!has_power(exponent)) {
            throw std::out_of_range(
                context + " requires unavailable P" + std::to_string(exponent)
            );
        }
    };

    const auto ensure_pure_routes = [&](const std::vector<RouteOp>& ops) {
        for (const RouteOp& op : ops) {
            if (op.kind != RouteKind::Pure)
                continue;

            const int target = op.lhs + op.rhs;
            if (has_power(target))
                continue;

            require_power(op.lhs, "pure route");
            require_power(op.rhs, "pure route");
            this->add_step(power_key(op.lhs), power_key(op.rhs), 'x');
        }
    };

    const auto serialize_monomial = [&] (
        double coefficient,
        int degree,
        bool multA,
        const std::vector<RouteOp>& ops
    ) -> std::string {
        ensure_pure_routes(ops);

        if (!multA) {
            require_power(degree, "pure monomial");
            if (coefficient == 1.0)
                return power_key(degree);

            const std::string coeff_key = this->add_value(
                'C',
                this->coeff_count++,
                std::format("{}", coefficient)
            );
            return this->add_step(coeff_key, power_key(degree), 'x');
        }

        const std::string coeff_key = this->add_value(
            'C',
            this->coeff_count++,
            std::format("{}", coefficient)
        );
        std::map<int, std::string> coefficient_terms;

        for (const RouteOp& op : ops) {
            if (op.kind == RouteKind::Pure)
                continue;

            if (op.kind == RouteKind::CoeffDirect) {
                require_power(op.rhs, "coeff_direct route");
                coefficient_terms[op.rhs] = this->add_step(
                    coeff_key,
                    power_key(op.rhs),
                    'x'
                );
                continue;
            }

            require_power(op.rhs, "coefficient route");

            std::string lhs_key;
            const auto lhs_it = coefficient_terms.find(op.lhs);
            if (lhs_it != coefficient_terms.end()) {
                lhs_key = lhs_it->second;
            } else {
                require_power(op.lhs, "coefficient base");
                lhs_key = this->add_step(
                    coeff_key,
                    power_key(op.lhs),
                    'x'
                );
                coefficient_terms[op.lhs] = lhs_key;
            }

            coefficient_terms[op.lhs + op.rhs] = this->add_step(
                lhs_key,
                power_key(op.rhs),
                'x'
            );
        }

        const auto found = coefficient_terms.find(degree);
        if (found != coefficient_terms.end())
            return found->second;

        require_power(degree, "coefficient-bearing monomial");
        return this->add_step(coeff_key, power_key(degree), 'x');
    };

    const auto serialize_terminal = [&]() -> std::string {
        std::string result_key;

        if (dcmp->term_plans.empty()) {
            for (int degree = 0;
                 degree < static_cast<int>(dcmp->coeff.size());
                 ++degree) {
                const double coefficient = dcmp->coeff[degree];
                if (coefficient == 0.0)
                    continue;

                std::string term_key;
                if (degree == 0) {
                    term_key = this->add_value(
                        'C',
                        this->coeff_count++,
                        std::format("{}", coefficient)
                    );
                } else {
                    term_key = serialize_monomial(
                        coefficient,
                        degree,
                        false,
                        {}
                    );
                }

                result_key = result_key.empty()
                    ? term_key
                    : this->add_step(result_key, term_key, '+');
            }
            return result_key;
        }

        for (const TermPlan& term : dcmp->term_plans) {
            if (term.degree < 0
                || static_cast<std::size_t>(term.degree) >= dcmp->coeff.size()) {
                throw std::out_of_range(
                    "term plan degree is outside polynomial coefficient range"
                );
            }

            const double coefficient = dcmp->coeff[term.degree];
            if (coefficient == 0.0)
                continue;

            std::string term_key;
            if (term.degree == 0) {
                term_key = this->add_value(
                    'C',
                    this->coeff_count++,
                    std::format("{}", coefficient)
                );
            } else {
                term_key = serialize_monomial(
                    coefficient,
                    term.degree,
                    term.multA,
                    term.route_ops
                );
            }

            result_key = result_key.empty()
                ? term_key
                : this->add_step(result_key, term_key, '+');
        }
        return result_key;
    };

    if (!dcmp->dcmp_p && !dcmp->dcmp_q)
        return serialize_terminal();

    // class_xi의 power는 자식 평가 전에 준비되어 있어야 한다.
    const double xi_coefficient = dcmp->xi.multA ? dcmp->xi.coeff : 1.0;
    const std::string xi_key = serialize_monomial(
        xi_coefficient,
        dcmp->xi.n,
        dcmp->xi.multA,
        dcmp->xi.route_ops
    );

    std::string p_key;
    std::string q_key;

    if (dcmp->eval_order == "q_then_p" || dcmp->eval_order == "q_only") {
        if (dcmp->dcmp_q)
            q_key = serialize_dcmp(dcmp->dcmp_q);
        if (dcmp->eval_order != "q_only" && dcmp->dcmp_p)
            p_key = serialize_dcmp(dcmp->dcmp_p);
    } else {
        if (dcmp->dcmp_p)
            p_key = serialize_dcmp(dcmp->dcmp_p);
        if (dcmp->eval_order != "p_only" && dcmp->dcmp_q)
            q_key = serialize_dcmp(dcmp->dcmp_q);
    }

    std::string result_key = xi_key;
    if (!p_key.empty())
        result_key = this->add_step(xi_key, p_key, 'x');
    if (!q_key.empty())
        result_key = this->add_step(result_key, q_key, '+');

    return result_key;
}

std::string EvalStep::get_value(std::string key)
{
    if (key.empty())
        throw std::invalid_argument("EvalStep::get_value: empty key");

    switch(key[0]) {
        case 'P': return powers.at(key);
        case 'C': return coeffs.at(key);
        case 'T': return terms.at(key);
        default:
            throw std::invalid_argument(
                "EvalStep::get_value: invalid key prefix: " + key
            );
    }
}

std::string EvalStep::add_value(
    char key_header,
    int key_number,
    std::string data
)
{
    const std::string key = std::string(1, key_header) + std::to_string(key_number);

    switch(key_header) {
        case 'C':
            this->coeffs.insert_or_assign(key, data);
            this->eval_step.push_back(step{'o', key, data, ""});
            break;
        case 'P':
            this->powers.insert_or_assign(key, data);
            if (key_number == 1)
                this->eval_step.push_back(step{'o', key, data, ""});
            break;
        case 'T':
            this->terms.insert_or_assign(key, data);
            break;
        default:
            throw std::invalid_argument("EvalStep::add_value: invalid key header");
    }
    return key;
}

std::string EvalStep::add_step(
    std::string key1,
    std::string key2,
    char op
)
{
    if (key1.empty() || key2.empty())
        throw std::invalid_argument("EvalStep::add_step: empty key");

    char key_header;
    int key_number;
    std::string data;

    if (op == '+') {
        key_header = 'T';
        key_number = this->term_count++;
        data = std::format("({} + {})", this->get_value(key1), this->get_value(key2));
    } else if (op == 'x') {
        if (key1[0] == 'P' && key2[0] == 'P') {
            const int lhs = std::stoi(key1.substr(1));
            const int rhs = std::stoi(key2.substr(1));
            key_header = 'P';
            key_number = lhs + rhs;
            data = std::format("x^{}", key_number);
        } else {
            const auto valid_prefix = [](char c) {
                return c == 'P' || c == 'C' || c == 'T';
            };
            if (!valid_prefix(key1[0]) || !valid_prefix(key2[0])) {
                throw std::invalid_argument(
                    std::format(
                        "EvalStep::add_step: invalid key data {} {} {}",
                        key1,
                        key2,
                        op
                    )
                );
            }

            key_header = 'T';
            key_number = this->term_count++;
            data = std::format("({} * {})", this->get_value(key1), this->get_value(key2));
        }
    } else {
        throw std::invalid_argument("EvalStep::add_step: invalid operation");
    }

    const std::string save_key = this->add_value(key_header, key_number, data);
    this->eval_step.push_back(step{op, key1, key2, save_key});
    return save_key;
}

void EvalStep::print_step()
{
    int count = 1;
    for (step s : this->eval_step) {
        this->print_line(s, count);
        ++count;
    }
}

void EvalStep::print_line(step s, int& count)
{
    std::string op_string;
    switch(s.op) {
        case '+': op_string = "ADD"; break;
        case 'x': op_string = (s.key1 == s.key2) ? "SQR" : "MUL"; break;
        case 'o': op_string = "RES"; break;
        default: op_string = "UNKNOWN"; break;
    }

    std::cout << std::format(
        "{}\t{}\t{}\t{}\t{}\n",
        count,
        op_string,
        s.key1,
        s.key2,
        s.save_key
    );
}
