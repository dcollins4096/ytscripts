
string = r'^Level\[(\d)\]: dt = ([^\s]+)\s*([^\s]+)\s*\(([^\s]+)/([^\s]+)\)'
s = 'Level[8]: dt = 6.3029e-05  0.000832312 (0.00172426/0.00172426)'
m = re.match(string, s)
print m
if m is not None:
    print m.groups()
