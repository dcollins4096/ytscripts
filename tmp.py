
for key in oned['full']:
    if key not in oneN['full']:
        print("Shit:",key)
    else:
        va = oned['full'][key]
        vb = oneN['full'][key]
    	print("ok", key,  va-vb)
for key in oneN['full']:
    if key not in oned['full']:
        print("Double Shit",key)

