execfile('go')
if 'yt' not in dir():
    import yt

ds = yt.load('/Users/dcollins/scratch/EnzoProjects/E12/jl01_JamesPancake1dUnigrid_CT/DD0000/jl010000')
ad = ds.all_data()
den = ad['density']
