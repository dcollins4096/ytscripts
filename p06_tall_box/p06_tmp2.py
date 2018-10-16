from go import *
def but1( data_1, data_2):
    plt.clf()
    p1=data_1['cooling_time'] > 0
    p2=data_2['cooling_time'] > 0
    m1=data_1['cooling_time'] < 0
    m2=data_2['cooling_time'] < 0
    plt.plot(data_1['temperature'][p1],data_1['cooling_time'][p1],c='r')
    plt.plot(data_1['temperature'][m1],-1*data_1['cooling_time'][m1],c='r', linestyle='--')
    plt.plot(data_2['temperature'][p2],data_2['cooling_time'][p2],c='k')
    plt.plot(data_2['temperature'][m2],-1*data_2['cooling_time'][m2],c='k', linestyle='--')
    plt.xscale('log'); plt.yscale('log')
    plt.ylabel('drooing_time '+str(data_2['cooling_time'].units))
    plt.savefig('p06_times.pdf')
