from IPython import embed
from collections import Counter

from sage.libs.singular.function_factory import ff

import resource
import sys
resource.setrlimit(resource.RLIMIT_STACK, [0x10000000, resource.RLIM_INFINITY])
sys.setrecursionlimit(0x100000)

attach('bapolar.py')
attach('misc.py')
attach('symmetry.py')
attach('hyperplane.py')
attach('tensors.py')

