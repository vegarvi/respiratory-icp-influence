from pylab import *
import scipy.interpolate as Spline
import scipy.signal as sc
from scipy.fftpack import fft


if __name__=='__main__':
	f1 = sys.argv[1]
	f2 = sys.argv[2]


def filters(filename1,filename2,cf=1./5,cfdiff=1./30,comma=True):
	if comma:
		p1 = []
		p2 = []
	
		f1 = open(filename1,'r')
		for line in f1:
			p1.append(float(line.replace(',','.')))
		f1.close()
		f2 = open(filename2,'r')
		for line in f2:
			a = float(line.replace(',','.'))
			p2.append(a)
	
		f2.close()
	
		p1 = array(p1)
		p2 = array(p2)
	else:
		p1 = loadtxt(filename1)
		p2 = loadtxt(filename2)
	#CONVENTION! P1 is inner pressure ICP CSF
	# P2 is outer pressure ICP EPI



	# Cutoff freq, normalized to Nyquist frequency, so cf=1 corresponds to 100 Hz. cf=1./15 used
	# USED 1./10 IN PAPER
	
	t = linspace(0,(len(p1)-1)*1./200,len(p1))
	b,a = sc.butter(7,cf)
	b2,a2 = sc.butter(7,cfdiff)
	dp = p2-p1
	dP = sc.filtfilt(b2,a2,p2-p1)
	P2 = sc.filtfilt(b,a,p2)
	P1 = sc.filtfilt(b,a,p1)
	return dp,p1,p2

'''
P1 = sc.filtfilt(b,a,p1)
P2 = sc.filtfilt(b,a,p2)
plot(t,P1)
show()
'''

'''
N = len(p1)
f = 1./200.
xf = linspace(0.0,1./(2.*f),N//2)
yf = fft(p1)
idx = where(xf<0.8)[0][-1]
idx2 = where(xf>1.4)[0][0]
print idx	
yf[0:idx]=0
yf[idx2:]=0
plot(yf)
show()
y2 = ifft(yf)
plot(abs(y2))
show()
'''
#plot(xf,2.0/N * np.abs(yf[0:N//2]))

#savetxt('p1_inner.txt',P1)
#savetxt('p2_outer.txt',P2)

if __name__=='__main__':
	dp,p1,p2 = filters(f1,f2)
	t = linspace(0,(len(dp)-1)*1./200,len(dp))
	N = len(p1)
	f = 1./200.
	xf = linspace(0.0,1./(2.*f),N//2)
	yf = fft(p1)
	plot(xf,abs(yf)[:N//2])
	show()
	'''
	if len(sys.argv[1])>2:
		savetxt('dP_P2_min_P1.txt',dp)
	plot(t,dp,linewidth=3)
	xlabel('time [s]')
	ylabel('dP [mmHg]')
	show()
	'''


	#plot(t,dp,t,p2-p1,'--')
	#show()
	
