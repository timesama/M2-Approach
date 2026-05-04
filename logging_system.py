import sys
from pathlib import Path


def _is_target_code(filename: str) -> bool:
    if not filename:
        return False
    p = Path(filename)
    parts = p.parts
    if "site-packages" in parts:
        return False
    # only app code we care about
    return any(token in parts for token in ("controllers", "dialogs")) or p.name in {"main.py", "main_window.py"}


def _format_file_lengths(frame):
    hints = []
    for name, value in frame.f_locals.items():
        lower = name.lower()
        if "file" not in lower:
            continue
        if isinstance(value, (list, tuple)):
            hints.append(f"{name}:len={len(value)}")
        elif isinstance(value, str):
            hints.append(f"{name}:str")
    return ", ".join(hints)


def _profile(frame, event, arg):
    if event != "call":
        return _profile
    code = frame.f_code
    filename = code.co_filename
    if not _is_target_code(filename):
        return _profile
    fn = code.co_name
    # ignore noisy internals / Qt event floods
    noisy_prefixes = ("event", "mouse", "paint", "resize", "enter", "leave", "timer")
    if fn.startswith("_") or fn.lower().startswith(noisy_prefixes):
        return _profile
    details = _format_file_lengths(frame)
    if details:
        print(f"[TRACE] {fn} | {details}")
    else:
        print(f"[TRACE] {fn}")
    return _profile


def install_function_logger():
    sys.setprofile(_profile)
