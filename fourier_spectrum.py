from pylab import *
import scipy.interpolate as Spline
import scipy.signal as sc
import sys
import os
from filtering import filters
from scipy.optimize import curve_fit
from convolution import *


def read(pat,freq=1./200,cut_hours=0):
	hours = cut_hours
	string = 'PatID%d'%pat
	for files in os.walk(string):    # if dp-file does not exist
		f = sorted(files[2])
		CSF = f[0]
		Ventricles = f[1]

	toss = 1./freq*2*60*60
	dp,p1,p2 = filters(string+'/'+CSF,string+'/'+Ventricles)
	if hours != 0:
		dp,p1,p2 = dp[toss:-toss*hours],p1[toss:-toss*hours],p2[toss:-toss*hours]
	daytime = CSF.split('_')[3]
	t0 = float(daytime[:2])*3600 + float(daytime[2:4])*60 + float(daytime[4:]) + 2*60*60
	t1 = t0 + (len(dp)-1)*freq
	T = (len(dp)-1)*freq
	t = linspace(0,T,len(dp))
	return daytime,t,dp,p1,p2

def fft_plot(dp,plt=True):
	dp-=average(dp)

	N = len(dp)
	f = 1./200
	xf = linspace(0,1./(2*f),N//2)
	y = fft(dp)
	yf = abs(y)[:N//2]*2.0/N
	idx = where(xf>2)[0][0]
	#plot(xf,yf/max(yf))
	#axis([0,1.7,0,0.5])	
	#xlabel('Frequency [Hz]')
	#ylabel('Amplitude (normalized)')
	b,a = sc.butter(7,0.15,btype='low')
	dP = sc.filtfilt(b,a,yf[:idx])
	if plt:
		#plot(xf[:idx],dP[:idx])
		show()
	return xf[:idx],dP[:idx]

def find_6mins(dp,minutes=6,threshold=1.5):
	dp -= average(dp)
	Smin = int(minutes*60*200)

	I = []
	idx = 0
	N = len(dp)
	
	while idx < N:
		start = idx
		cur = idx
		M0 = dp[idx]
		M1 = dp[idx]

		while abs(M1-M0)<threshold and idx<N:
			if dp[idx]>M1:
				M1 = dp[idx]
				#print 'Found Max: %g at index %g'%(M1,idx)
				#raw_input()
			if dp[idx]<M0:
				M0 = dp[idx]
				#print 'Found Min: %g at index %g'%(M0,idx)
				#raw_input()
			cur+=1
			idx+=1

			if cur-start >= Smin:
				#print 'success in [%g, %g]'%(start,cur-1)
				I.append([start,cur-200])
				break;

			
	return array(I)


def return_amps(dp,I,window_size=10,pl=True,denoise=True,resp=[0.2,0.4],card=[0.7,1.6]):
	w = window_size
	I = array(I)
	I_out = []
	Ramp = []
	Camp = []
	Rfreq = []
	Cfreq = []
	for i in range(len(I)):
		'''
		print I[i][0]
		print I[i][1]
		print I
		'''
		func = dp[I[i][0]:I[i][1]]
		#plot(func)
		#show()
		
		xf,fftp = fft_plot(func,plt=False)
		idx = where(xf>0.6)[0][0]
		if denoise:
			popt, pcov = curve_fit(noise,xf[:idx],fftp[:idx])
			if pl:
				#subplot(2,1,1)
				#plot(xf,fftp,xf,noise(xf,*popt))
				#subplot(2,1,2)
				plot(xf,fftp-noise(xf,*popt),'k',linewidth=2)
				xlabel('Frequency [Hz]')
				ylabel('Amplitude')
				legend(['fft(dICP)'],loc='upper left')
				show()
			fftp -= noise(xf,*popt)
		idx0 = where(xf>resp[0])[0][0]
		idx1 = where(xf>resp[1])[0][0]
		idx2 = where(xf>card[0])[0][0]
		idx3 = where(xf>card[1])[0][0]
		#print idx0,idx1,idx2,idx3
		Ridx = where(fftp[idx0:idx1] == max(fftp[idx0:idx1]))[0][0]
		Cidx = where(fftp[idx2:idx3] == max(fftp[idx2:idx3]))[0][0]
		if w==0:
			mr = fftp[idx0+Ridx]
			mc = fftp[idx2+Cidx]
		else:		
			mr = average(fftp[idx0+Ridx-w:idx0+Ridx+w])
			mc = average(fftp[idx2+Cidx-w:idx2+Cidx+w])
		if mc > 1e-3:
			Ramp.append(mr)
			Camp.append(mc)
			Rfreq.append(xf[idx0+Ridx])
			Cfreq.append(xf[idx2+Cidx])
			I_out.append(I[i])
			#print mc, mr
	#print 'N = ', len(Ramp)
	
	return array(Ramp),array(Camp),array(Rfreq),array(Cfreq), I_out

def noise(t,A0,tau):
	return A0*exp(-t/tau)

if __name__ =='__main__':
	matplotlib.rc('font',size=20)


	daytime,t,dp,_,_ = read(int(sys.argv[1]),cut_hours=0)
	t0 = float(daytime[:2])*3600 + float(daytime[2:4])*60 + float(daytime[4:]) + 2*60*60
	sleep = 24*60*60 + 1*60*60 - t0 # from 1 a.m.
	wake = 24*60*60 + 5*60*60  - t0 # to 5 a.m.
	minutes = 6




	#p = dp[:1e7]
	#plot(dp)
	#show()
	treshold = 2.
	I0 = find_6mins(dp,minutes,treshold)
	print(('Recording startime = ',daytime)) # %d:%d:%d'%(t0//3600,t0%3600//60,(t0%3600)%60)
	
	
	#tmax,pmax,tmin,pmin = convolution(dp)

	R,C,Rf,Cf,I = return_amps(dp,array(I0),window_size=0,pl=True)
	T = []
	for i in range(len(I)):
		T.append(t[I[i][0]])
	print( T)
	print(( 'Respiration = %.3e +- %.3e     f = %.3f'%(average(R),std(R,ddof=1),average(Rf))))
	print(( 'Cardiac = %.3e +- %.3e     f = %.3f'%(average(C),std(C,ddof=1),average(Cf))))

	print(( '\n RATIO: ',average(C)/average(R)))
	print(( len(T),len(Rf)))
	plot(T,Rf,'ro')
	plot(T,Cf,'bo')
	plot([sleep,sleep+1e-10],[-100,100],'--k')
	plot([wake,wake+1e-10],[-100,100],'--k')
	xlabel('Time [s]')
	ylabel('ICP difference')
	axis([T[0],T[-1],min(Rf),max(Cf)])
	legend(['Respiratory', 'Cardiac'])
	show()


# Very good examples of respiration: 
# Pat1 2 880 000 - 3 020 000


# Possibly:
# Pat3 3025000 - 3050000
# Pat4 9111000 - 9127000
# Pat8 

# CHiari Pat:
# Pat12172013 58000 - 70000
