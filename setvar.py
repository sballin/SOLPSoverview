#!/bin/env python3
import sys
import re

def main(filename: str, var_name: str, new_value: str):
    variable_pattern = re.compile(rf"'{var_name}'\s*'([^']*)'")
    with open(filename, 'rb') as f:
        lines = f.readlines()
    for i, line in enumerate(lines):
        line = line.decode('utf8')
        match = variable_pattern.search(line)
        if match:
            start, end = match.span(1)  # find the span of group 1 which is the value
            lines[i] = (line[:start] + new_value + line[end:]).encode('utf8')  # replacements
    with open(filename, 'wb') as f:  # overwrite the file 
        f.writelines(lines)


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print('Usage: setvar ./b2mn.dat b2mndr_dtim 1e-3')
        sys.exit()
    filename = sys.argv[1]
    var_name = sys.argv[2]
    new_value = sys.argv[3]
    main(filename, var_name, new_value)