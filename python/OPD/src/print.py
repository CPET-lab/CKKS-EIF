# print.py

def poly_to_str(poly: list[float]) -> str:
    """
    Converts polynomial into string.
    """
    res = ""
    space = ""
    for i in range(len(poly) - 1, -1, -1):
        coeff = poly[i]
        if coeff == 0:
            continue

        # sign
        if not res:
            mark = "-" if coeff < 0 else ""
        else:
            mark = "-" if coeff < 0 else "+"

        c = abs(coeff)

        # numeric formatting
        if float(c).is_integer():
            c_str = str(int(c))
        else:
            c_str = repr(c)

        # term formatting (omit 1 for non-constant terms)
        if i == 0:
            term = c_str
        elif i == 1:
            term = ("x" if c == 1 else f"{c_str}x")
        else:
            term = (f"x^{i}" if c == 1 else f"{c_str}x^{i}")

        res += f"{space}{mark}{term}"
        space = " "

    return res or "0"


def type_decomp_to_str(dcmp) -> str:
    """Print a type decomposition using original coefficient symbols.

    Each child stores its own representative coefficients, so numerical child
    values must not be used to display an identity for the original polynomial.
    """
    coefficients = [
        "0" if kind == "0" else f"a{degree}"
        for degree, kind in enumerate(dcmp.coeff_type)
    ]

    def power(degree):
        return "x" if degree == 1 else f"x^{degree}"

    def polynomial(values):
        terms = []
        for degree in range(len(values) - 1, -1, -1):
            coefficient = values[degree]
            if coefficient == "0":
                continue
            if degree == 0:
                terms.append(coefficient)
            elif coefficient == "1":
                terms.append(power(degree))
            else:
                terms.append(f"({coefficient})*{power(degree)}")
        return " + ".join(terms) or "0"

    def render(node, values):
        if node.xi.n == 0:
            return polynomial(values)

        degree = node.xi.n
        quotient = values[degree:]
        factor = power(degree)
        if node.xi.multA:
            leading = values[-1]
            factor = f"({leading})*{factor}"
            quotient = [
                "0" if value == "0" else "1" if value == leading
                else value if leading == "1" else f"({value})/({leading})"
                for value in quotient
            ]
        if node.dcmp_p is not None:
            factor = f"{factor}*({render(node.dcmp_p, quotient)})"
        if node.dcmp_q is not None:
            factor += f" + ({render(node.dcmp_q, values[:degree])})"
        return factor

    return render(dcmp, coefficients)
