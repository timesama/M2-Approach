import os
import sys
from pathlib import Path


def _is_repo_code(filename: str) -> bool:
    if not filename:
        return False
    parts = Path(filename).parts
    skip = {"site-packages", "lib-dynload", "importlib", "encodings"}
    return not any(p in skip for p in parts)


def _format_file_lengths(frame):
    hints = []
    for name, value in frame.f_locals.items():
        lower = name.lower()
        if "file" not in lower:
            continue
        if isinstance(value, (list, tuple)):
            hints.append(f"{name}:len={len(value)}")
        elif isinstance(value, str):
            size = None
            if os.path.exists(value) and os.path.isfile(value):
                try:
                    size = os.path.getsize(value)
                except OSError:
                    size = None
            if size is not None:
                hints.append(f"{name}:path,size={size}")
            else:
                hints.append(f"{name}:str")
    return ", ".join(hints)


def _profile(frame, event, arg):
    if event != "call":
        return _profile
    code = frame.f_code
    filename = code.co_filename
    if not _is_repo_code(filename):
        return _profile
    fn = code.co_name
    details = _format_file_lengths(frame)
    if details:
        print(f"[TRACE] {fn} | {details}")
    else:
        print(f"[TRACE] {fn}")
    return _profile


def install_function_logger():
    sys.setprofile(_profile)
