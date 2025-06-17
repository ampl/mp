import os
import re
import sys

# Paths to the relevant files (relative to this script)
script_dir = os.path.dirname(os.path.abspath(__file__))
root = os.path.join(script_dir, "..")
base_dir = os.path.join(root, "solvers")


def get_latest_changelog_entry(changelog_path):
    """Extracts and returns the latest changelog entry (date and content) from a markdown changelog file."""
    with open(changelog_path, encoding="utf-8") as f:
        lines = f.readlines()
    entry_lines = []
    in_entry = False
    for i, line in enumerate(lines):
        if line.startswith("## "):
            if in_entry:
                break  # End of the latest entry
            entry_lines = [line.rstrip()]
            in_entry = True
        elif in_entry:
            entry_lines.append(line.rstrip())
    # Extract date
    date = None
    for l in entry_lines:
        m = re.match(r"^## (\d+)", l.strip())
        if m:
            date = int(m.group(1))
            break
    return entry_lines, date


def propagate_mp_update(names, mp_entry, mp_date):
    """For each name, update CHANGES.name.md in base_dir/name if its latest entry is older than mp_entry."""
    formatted_mp_entry = format_mp_entry_with_bullet(mp_entry)
    print("\nThe following entry will be inserted:\n")
    print(formatted_mp_entry)
    confirm = input("\nProceed with updating all changelogs? [y/N]: ").strip().lower()
    if confirm != 'y':
        print("Aborted by user.")
        return
    for name in names:
        changes_path = os.path.join(base_dir, name, f"CHANGES.{name}.md")
        if not os.path.exists(changes_path):
            print(f"File not found: {changes_path}")
            continue
        latest_entry, latest_date = get_latest_changelog_entry(changes_path)
        if latest_date is not None and latest_date >= mp_date:
            print(f"{name}: up to date (latest {latest_date} >= mp {mp_date})")
            continue
        # Insert formatted mp_entry after the header
        with open(changes_path, encoding="utf-8") as f:
            lines = f.readlines()
        header_end = 0
        for i, line in enumerate(lines):
            if line.strip() == "" and i > 0:
                header_end = i
                break
        new_lines = lines[:header_end] + [formatted_mp_entry + "\n\n"] + lines[header_end:]
        with open(changes_path, "w", encoding="utf-8") as f:
            f.writelines(l + ("\n" if not l.endswith("\n") else "") for l in new_lines)
        print(f"{name}: updated with mp entry {mp_date}")


def format_mp_entry_with_bullet(mp_entry):
    header = mp_entry[0]
    # Everything after the header is the content
    content = mp_entry[1:]
    # Remove leading/trailing blank lines in content
    while content and not content[0].strip():
        content = content[1:]
    while content and not content[-1].strip():
        content = content[:-1]
    # Indent all content lines by 2 spaces
    indented = ["  " + l if l.strip() else "" for l in content]
    # Add the bullet point
    result = ["", "", header, "- Changes in MP"] + indented
    return "\n".join(result)


opensource = ["baronmp", "cbcmp", "gcgmp", "scipmp"]
mip = ["copt", "cplex", "gurobi", "highsmp", "mosek", "xpress"]


if __name__ == "__main__":
    changelog_path = os.path.join(script_dir, "..", "CHANGES.mp.md")
    latest_entry, latest_date = get_latest_changelog_entry(changelog_path)
    propagate_mp_update(opensource +mip, latest_entry, latest_date)

