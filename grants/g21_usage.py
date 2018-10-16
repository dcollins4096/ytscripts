if 'ef' not in dir():
    execfile('go')
import pdb
fname = '/Users/dcollins/RESEARCH2/Paper49_EBQU/2017-09-14/_Tracking_the_Assembly_History_of_Star_Formation_through_Clouds_2017-09-15T15-02-41.530Z.csv'
ef=execfile
lines = open(fname,'r').readlines()
# "Project","Resource","User","Start Time","End Time", "SUs" 
def date_scrub(date):
    full_day, time , timezone = date.split(" ")
    month, day,year = map(int,full_day.split("/"))
    hour, minute = map(int,time.split(":"))
    dt = datetime.datetime(year=year,month=month,day=day,hour=hour,minute=minute)
    return dt

class job():
    def __init__(self,line):
        split = [f.strip('"') for f in line.split(",")]
        if len(split) == 7:
            self.Project = split[0]
            self.Resource = split[1]
            self.User = split[2]
            self.Start = date_scrub(split[3])
            self.End = date_scrub(split[4])
            self.SUs = float(split[5])
            self.ok = True
        else:
            print "Wrong line length.", len(split)
            self.ok = False
job_array = []
bad_line_counter = 0
for line in lines[1:]:
    this_job = job(line)
    if this_job.ok:
        job_array.append( job(line) )
    else:
        bad_line_counter += 1

#print "bad lines:", bad_line_counter
SUs = [j.SUs for j in job_array]
starts = [j.Start for j in job_array]
#print min(SUs), max(SUs)

plt.clf()
plt.hist(SUs,histtype='step')
plt.savefig('g21_xsede_sus.pdf')
plt.clf()
plt.plot(starts,np.cumsum(SUs))
plt.savefig('g21_cumulative.pdf')
