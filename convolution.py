from filtering import filters
from pylab import *
import scipy.signal as sc

if __name__=='__main__':

	P = loadtxt(sys.argv[1])

def convolution(P,dist=80,minBCL=0.5,):
	Pmax = [0]
	tmax = [0]
	Pmin = [0]
	tmin = [0]
	P2max = [0]
	t2max = [0]
	P2min = [0]
	t2min = [0]
	Hz = 200.
	t = linspace(0,(len(P)-1)/Hz,len(P))


	FoundPeak = False

	for i in range(dist,len(P)-dist):
		if P[i] >= max([P[j] for j in range(i-dist,i)]) and P[i] >= max([P[j] for j in range(i+1,i+dist)]):
			if abs(tmax[-1]-t[i]) > minBCL:
				if FoundPeak==False:
					Pmax.append(P[i])
					tmax.append(t[i])
					FoundPeak = True
		if P[i] <= min([P[j] for j in range(i-dist,i)]) and P[i] <= min([P[j] for j in range(i+1,i+dist)]):
			if abs(tmin[-1]-t[i]) > minBCL:
				if FoundPeak==True:
					Pmin.append(P[i])
					tmin.append(t[i])
					FoundPeak=False

	return array(tmax),array(Pmax),array(tmin),array(Pmin)

if __name__=='__main__':
	cf = 1/30.
	b,a = sc.butter(7,cf)
	P = sc.filtfilt(b,a,P)
	tmax,Pmax,tmin,Pmin = convolution(P)
	t = linspace(0,(len(P)-1)/200.,len(P))
	plot(t,P)
	plot(tmin,Pmin,'o',tmax,Pmax,'o')
	print(average(Pmax)-average(Pmin))
	print(average(tmax)-average(tmin))
	Nmax = len(Pmax)
	Nmin = len(Pmin)
	if Nmax > Nmin: 
		N = Nmin
		dP = Pmax[1:] - Pmin
	elif Nmin > Nmax:
		dP = Pmax - Pmin[:-1]
		N = Nmax
	else:
		dP = Pmax - Pmin
	N = len(dP)
	DP = []
	T = []
	for i in range(N):
		if dP[i]<10:
			T.append(t[i])
			DP.append(dP[i])
		else:
			pass
	print(average(array(DP)))
	print(len(DP), len(dP))
	show()

