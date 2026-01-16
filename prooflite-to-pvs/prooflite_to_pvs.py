#!/usr/bin/env python3
"""
Convert ProofLite proof scripts to raw PVS commands.

This tool performs the following transformations:
1. Removes all `(THEN ` prefixes
2. Removes all `(SPREAD ` prefixes
3. Removes extraneous closing parentheses that become unbalanced
4. Strips the `%|- ` comment prefix from lines
"""

import re
import sys
from pathlib import Path


def strip_prooflite_prefix(text: str) -> str:
    """Remove the %|- prefix from ProofLite lines."""
    lines = text.split('\n')
    stripped_lines = []
    for line in lines:
        # Remove %|- prefix (with optional leading whitespace)
        stripped = re.sub(r'^(\s*)%\|-\s?', r'\1', line)
        stripped_lines.append(stripped)
    return '\n'.join(stripped_lines)


def remove_then_spread(text: str) -> str:
    """Remove (THEN and (SPREAD keywords."""
    # Remove (THEN with optional whitespace after
    text = re.sub(r'\(THEN\s+', '', text)
    # Remove (SPREAD with optional whitespace after
    text = re.sub(r'\(SPREAD\s+', '', text)
    return text


def balance_parentheses(text: str) -> str:
    """
    Remove extraneous closing parentheses by parsing the structure.

    After removing THEN and SPREAD, there will be unmatched closing parens.
    We need to remove them while preserving the structure.
    """
    result = []
    paren_depth = 0
    i = 0
    in_string = False

    while i < len(text):
        char = text[i]

        # Handle string literals (to avoid counting parens inside strings)
        if char == '"' and (i == 0 or text[i-1] != '\\'):
            in_string = not in_string
            result.append(char)
            i += 1
            continue

        if in_string:
            result.append(char)
            i += 1
            continue

        if char == '(':
            paren_depth += 1
            result.append(char)
        elif char == ')':
            if paren_depth > 0:
                paren_depth -= 1
                result.append(char)
            # else: skip the extraneous closing paren
        else:
            result.append(char)

        i += 1

    return ''.join(result)


def clean_whitespace(text: str) -> str:
    """Clean up excessive whitespace while preserving structure."""
    # Remove trailing whitespace on each line
    lines = text.split('\n')
    lines = [line.rstrip() for line in lines]

    # Remove empty lines at the start and end
    while lines and not lines[0].strip():
        lines.pop(0)
    while lines and not lines[-1].strip():
        lines.pop()

    return '\n'.join(lines)


def convert_prooflite_to_pvs(prooflite_text: str) -> str:
    """
    Convert a ProofLite proof script to raw PVS commands.

    Args:
        prooflite_text: The ProofLite proof script text

    Returns:
        The converted PVS commands
    """
    # Step 1: Strip the %|- prefix
    text = strip_prooflite_prefix(prooflite_text)

    # Step 2: Remove (THEN and (SPREAD
    text = remove_then_spread(text)

    # Step 3: Balance parentheses by removing extraneous closing parens
    text = balance_parentheses(text)

    # Step 4: Clean up whitespace
    text = clean_whitespace(text)

    return text


def main():
    """Main entry point for the CLI."""
    if len(sys.argv) < 2:
        print("Usage: python prooflite_to_pvs.py <input_file> [output_file]")
        print("       python prooflite_to_pvs.py -  (read from stdin)")
        sys.exit(1)

    input_arg = sys.argv[1]

    # Read input
    if input_arg == '-':
        prooflite_text = sys.stdin.read()
    else:
        input_path = Path(input_arg)
        if not input_path.exists():
            print(f"Error: Input file '{input_arg}' not found")
            sys.exit(1)
        prooflite_text = input_path.read_text()

    # Convert
    pvs_commands = convert_prooflite_to_pvs(prooflite_text)

    # Write output
    if len(sys.argv) >= 3:
        output_path = Path(sys.argv[2])
        output_path.write_text(pvs_commands)
        print(f"Output written to {output_path}")
    else:
        print(pvs_commands)


if __name__ == '__main__':
    main()
