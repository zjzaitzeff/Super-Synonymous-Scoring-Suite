import sys
import os
output1 = sys.argv[1]
output2 = sys.argv[2]
output3 = sys.argv[3]


if __name__ == '__main__':
    lst = [output1, output2, output3]
    for i in lst:
        if not os.path.exists(i):
            # Create parent directories if needed
            os.makedirs(os.path.dirname(i), exist_ok=True)

        # Create the empty file
            with open(i, 'w') as f:
                pass
            print(f"Created: {i}")
        else:
            print(f"Exists: {i}")
    