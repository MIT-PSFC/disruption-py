import cProfile
import pstats 
from pstats import SortKey

from disruption_py.shots import D3DShot

cProfile.run('D3DShot(175552, "EFIT05")', 'd3d_shot_stats')
p = pstats.Stats('d3d_shot_stats')
p.sort_stats(SortKey.CUMULATIVE).print_stats(20)
p.sort_stats(SortKey.TIME).print_stats(20)

