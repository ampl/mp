#!/usr/bin/env python3
"""Lay out AMPL models to be read: at most 80 characters per line, continuation
lines indented, blank lines between declarations, objectives and constraints.

Changes layout only, plus one token cleanup: a parenthesized numeric
coefficient `(c)*` loses its parentheses, `+ (-c)*` becomes `- c*`, and a
coefficient 1 is dropped. Statements longer than WIDTH are broken before
the operator of lowest precedence at the top level (logical, then
relational, then + and -), recursively; a long call or parenthesized group
puts one argument per line. Blank lines separate the declarations, the
objectives and the constraints.

Usage:
  amplfmt.py FILE.mod...            print the formatted models
  amplfmt.py --in-place FILE.mod... rewrite them
Check with nlcheck.py that the NL files do not change.
"""
import re
import sys

WIDTH = 80
INDENT = 3
KEYWORDS = {"if", "then", "else", "and", "or", "not", "in", "sum", "numberof",
            "count", "atleast", "atmost", "exactly", "alldiff", "exists",
            "forall", "min", "max", "abs", "prod"}
# break classes, lowest precedence first
CLASSES = [{"<==>", "==>", "else", "then"}, {"or", "||"}, {"and", "&&"},
           {"<=", ">=", "=", "==", "!=", "<", ">", "<>"}, {"+", "-"}]
TOKEN = re.compile(r"'[^']*'|\"[^\"]*\"|<==>|==>|:=|\*\*|<=|>=|!=|==|<>|&&|\|\||\.\.|"
                   r"[A-Za-z_][A-Za-z_0-9]*|(?:\d+(?:\.\d*)?|\.\d+)(?:[eEdD][-+]?\d+)?|\S")
NUM = r"\d+(?:\.\d+)?(?:[eE][-+]?\d+)?"
OPEN, CLOSE = "([{", ")]}"
AT_SPACE_ONLY = False      # careful mode: no '(' / ',' breaks; fill only at space


def clean_coefficients(s):
    """Only a coefficient that starts a term: after the start, ':', '(',
    ',' or a binary + or - (never after '^', '*' or '/')."""
    start = r"(^|[:(,]\s*|[^-+*/^(,:\s]\s*[+\-]\s*)"
    s = re.sub(r"\+\s*\(-(%s)\)\s*\*" % NUM, r"- \1*", s)
    s = re.sub(start + r"\((-?%s)\)\s*\*" % NUM, r"\1\2*", s)
    s = re.sub(start + r"(-?)1\*(?=[A-Za-z(])", r"\1\2", s)
    return s


def tokenize(s):
    out, pos = [], 0
    for m in TOKEN.finditer(s):
        out.append((m.group(), m.start() > pos))
        pos = m.end()
    return out


def flat(toks):
    return "".join((" " if sp and i else "") + t for i, (t, sp) in enumerate(toks))


def depths(toks):
    d, out = 0, []
    for t, _ in toks:
        if t in CLOSE:
            d -= 1
        out.append(d)
        if t in OPEN:
            d += 1
    return out


def is_binary(toks, i):
    if i == 0:
        return False
    p = toks[i - 1][0]
    return p in CLOSE or (p[0].isalnum() or p[0] == "_") and p not in KEYWORDS


def split_at(toks, pos):
    segs, last = [], 0
    for p in pos:
        segs.append(toks[last:p])
        last = p
    segs.append(toks[last:])
    return [s for s in segs if s]


def layout(toks, ind):
    s = flat(toks)
    if ind + len(s) <= WIDTH:
        return [" " * ind + s]
    d = depths(toks)
    for cls in CLASSES:
        pos = [i for i, (t, sp) in enumerate(toks) if i and d[i] == 0 and t in cls
               and (t not in "+-" or is_binary(toks, i))
               ]   # (careful mode too: a break before an operator adds only whitespace)
        if pos:
            segs = split_at(toks, pos)
            if "+" in cls:              # sums: several short terms per line
                segs = pack(segs, ind)
            return [ln for seg in segs for ln in layout(seg, ind)]
    if not AT_SPACE_ONLY:
        return layout_group(toks, ind, d) or [" " * ind + s]
    return paren_group(toks, ind, d) or fill(toks, ind)


def paren_group(toks, ind, d):
    """Careful mode: 'prefix(inner)suffix' with a short prefix and suffix
    (a sign, 'sum {...}', a function name): the inner part is laid out
    under the '(' (no break is added next to the parentheses)."""
    for k, (t, _) in enumerate(toks):
        if t != "(" or d[k] != 0:
            continue
        m = matching(toks, k, d)
        head = flat(toks[:k + 1])
        if m is None or m < k + 2 or len(head) > 40 or len(flat(toks[m:])) > 20:
            continue
        tail = flat([(toks[m][0], False)] + toks[m + 1:])
        global WIDTH
        WIDTH, width = WIDTH - len(tail), WIDTH     # the last line keeps room
        try:
            inner = layout(toks[k + 1:m], ind + len(head))
        finally:
            WIDTH = width
        inner[0] = " " * ind + head + inner[0].lstrip()
        inner[-1] += tail
        return inner
    return None


def fill(toks, ind):
    """Greedy fill, breaking before tokens preceded by space or before a
    binary operator at any depth (only whitespace is added). On overflow
    the line is broken at its last break point."""
    ops = set().union(*CLASSES) | {"*", "/"}
    lines, cur, last = [], [], None        # last: break point within cur
    for i, t in enumerate(toks):
        if cur and (t[1] or (t[0] in ops and is_binary(toks, i))):
            last = len(cur)
        cur.append(t)
        if ind + len(flat(cur)) > WIDTH and last:
            lines.append(" " * ind + flat(cur[:last]))
            cur, last = cur[last:], None
    return lines + [" " * ind + flat(cur)]


def pack(segs, ind):
    out = []
    for sg in segs:
        if out and ind + len(flat(out[-1] + sg)) <= WIDTH:
            out[-1] = out[-1] + sg
        else:
            out.append(sg)
    return out


def matching(toks, k, d):
    for j in range(k + 1, len(toks)):
        if toks[j][0] in CLOSE and d[j] == d[k]:
            return j
    return None


def layout_group(toks, ind, d):
    """A parenthesized group or a call, the one that ends last at depth 0."""
    ks = [k for k, (t, _) in enumerate(toks) if t == "(" and d[k] == 0]
    if not ks:
        return None
    k = ks[-1] if toks[-1][0] != ")" else ks[0]
    m = matching(toks, k, d)
    if m is None:
        return None
    if k == 0 and m == len(toks) - 1:            # (expr): align under "("
        inner = layout(toks[1:-1], ind + 1)
        inner[0] = " " * ind + "(" + inner[0].lstrip()
        inner[-1] += ")"
        return inner
    head, body, tail = toks[:k + 1], toks[k + 1:m], toks[m:]
    bd = depths(body)
    commas = [i for i, (t, _) in enumerate(body) if t == "," and bd[i] == 0]
    lines = [" " * ind + flat(head)]
    last = 0
    for c in commas + [len(body)]:
        arg = body[last:c]
        last = c + 1
        sub = layout(arg, ind + INDENT)
        if c < len(body):
            sub[-1] += ","
        lines += sub
    lines[-1] += flat([(tail[0][0], False)] + tail[1:])
    return lines


def kind(stmt):
    w = stmt.split(None, 1)[0] if stmt.strip() else ""
    if w in ("minimize", "maximize"):
        return "obj"
    if w in ("s.t.", "subject") or re.match(r"[A-Za-z_]\w*\s*(\{[^}]*\})?\s*:", stmt):
        return "con"
    return w


def format_statement(stmt):
    stmt = " ".join(stmt.split())
    head, sep, body = stmt.partition(":")
    if not sep or kind(stmt) not in ("obj", "con"):
        return [stmt] if len(stmt) <= WIDTH else layout(tokenize(stmt), 0)
    body = clean_coefficients(body.strip().rstrip(";").strip())
    one = f"{head}: {body};"
    if len(one) <= WIDTH:
        return [one]
    # ';' is a token of the body, so that the last line keeps room for it
    return [head + ":"] + layout(tokenize(body + ";"), INDENT)


def statements(text):
    """Comment lines and statements (to ';' outside brackets) in order."""
    cur, depth = "", 0
    for line in text.splitlines():
        if not cur.strip() and line.lstrip().startswith("#"):
            yield line.rstrip()
            continue
        for ch in line:
            cur += ch
            depth += (ch in OPEN) - (ch in CLOSE)
            if ch == ";" and depth == 0:
                yield cur.strip()
                cur = ""
        cur += " "
    if cur.strip():
        yield cur.strip()


def wrap_comment(line):
    if len(line) <= WIDTH:
        return [line]
    pre = line[:len(line) - len(line.lstrip())] + "# "
    words, out = line.lstrip()[1:].split(), [pre]
    for w in words:
        if len(out[-1]) + len(w) > WIDTH and out[-1] != pre:
            out.append(pre)
        out[-1] += ("" if out[-1] == pre else " ") + w
    return out


def format_model(text):
    out, prev = [], None
    for st in statements(text):
        if not st.strip():
            continue                             # blank lines are re-made
        if st.lstrip().startswith("#"):
            out += wrap_comment(st)
            continue
        k = kind(st)
        if prev is not None and k != prev and out and out[-1] != "" \
                and not out[-1].lstrip().startswith("#"):
            out.append("")
        out += format_statement(st)
        prev = k
    return "\n".join(out) + "\n"


def code_and_comment(line):
    """Split a line at a '#' outside quotes."""
    q = None
    for i, ch in enumerate(line):
        if q:
            q = None if ch == q else q
        elif ch in "'\"":
            q = ch
        elif ch == "#":
            return line[:i], line[i:]
    return line, ""


def wrap_line(line):
    """Careful mode for one physical line: continuation lines get the
    line's own indentation plus INDENT; a line with a comment, or a
    comment line, is left as it is."""
    global WIDTH
    code, comment = code_and_comment(line)
    if len(line.rstrip().expandtabs(8)) <= WIDTH or comment or not code.strip():
        return [line]
    lead = code[:len(code) - len(code.lstrip())]
    col = len(lead.expandtabs(8))
    WIDTH, width = WIDTH - INDENT, WIDTH      # room for the continuation indent
    try:
        out = layout(tokenize(code.strip()), col)
    finally:
        WIDTH = width
    res = []
    for i, ln in enumerate(out):
        pos = len(ln) - len(ln.lstrip())
        # continuations at the line's own column get INDENT; deeper lines
        # (aligned under a parenthesis) keep their column
        pad = pos - col + (INDENT if i and pos == col else 0)
        res.append(lead + " " * pad + ln.lstrip())
    return res


def wrap_lines(text):
    """Careful mode for third-party and converted models: only long code
    lines are broken, before operators or where they have whitespace;
    nothing but whitespace changes."""
    global AT_SPACE_ONLY
    AT_SPACE_ONLY, old = True, AT_SPACE_ONLY
    try:
        lines = [w for ln in text.split("\n") for w in wrap_line(ln)]
    finally:
        AT_SPACE_ONLY = old
    return "\n".join(lines)


def main(argv):
    inplace = argv[:1] == ["--in-place"]
    for f in argv[1:] if inplace else argv:
        text = open(f).read()
        new = format_model(text)
        if inplace:
            if new != text:
                open(f, "w").write(new)
        else:
            sys.stdout.write(new)


if __name__ == "__main__":
    main(sys.argv[1:])
