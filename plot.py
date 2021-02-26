#%%
import numpy as np 
import matplotlib.pyplot as plt
from math import floor
#%%
def load_dataset(filename):
    datas=np.loadtxt(filename)
    datas=datas.T
    print(datas.shape)
    datareal=datas[0::2,:]
    dataimage=datas[1::2,:] 
    datareal=datareal.astype(np.complex)
    datareal[1:,:]=datareal[1:,:]+1j*dataimage[1:,:]
    return datareal 

#%% 
file1="rf_specs.txt"
rf_specs=load_dataset(file1)
file2="corr_specs.txt"
corr_specs=load_dataset(file2)

#%%
shape=rf_specs.shape
leng=shape[0]
dt=0.01
df=1/dt/leng/2
f1=0
f2=3
ii1=floor(f1/df)
ii2=floor(f2/df)
freqs=np.linspace(0,(leng-1)*df,leng)
#%%
rf_specs_amp=abs(rf_specs)
rf_specs_amp_ave=rf_specs_amp[:,].mean(axis=1)
corr_specs_amp=abs(corr_specs)
corr_specs_amp_ave=corr_specs_amp[:,].mean(axis=1)
#%%
nline=shape[1]
fig, axs=plt.subplots(nrows=2,ncols=1,figsize=(6,8))
for i in range(2):

    axs[0].plot(freqs[ii1:ii2],corr_specs_amp[ii1:ii2,i],lw=0.3,c='g')
    axs[1].plot(freqs[ii1:ii2],rf_specs_amp[ii1:ii2,i],lw=0.3,c='b')
#set x label
axs[0].plot(freqs[ii1:ii2],corr_specs_amp_ave[ii1:ii2],lw=0.5,c='r')
axs[1].plot(freqs[ii1:ii2],rf_specs_amp_ave[ii1:ii2],lw=0.5,c='r')


axs[0].set_xlabel("freqs(hz)")
axs[1].set_xlabel("freqs(hz)")
#set y label
axs[0].set_ylabel("amp")
axs[1].set_ylabel("amp")
#set title
axs[0].set_title("amplitute of corr")
axs[1].set_title("amplitute of rf")
axs[1].set_ylim(0,1)
fig.tight_layout()
plt.show()

# %%
from scipy.fftpack import ifft
rf_specs_ave=rf_specs[:,70:90].mean(axis=1)
fc=2
imax=floor(fc/df)
rf_specs_ave[imax:]=0
hanning=np.linspace(0,imax-1,imax)
hanning=np.cos(np.pi*hanning/(imax-1)/2)**2
rf_specs_ave[:imax]=rf_specs_ave[:imax]*hanning
rf_specs_ave_conj=np.flip(np.conjugate(rf_specs_ave[1:]))
rf_specs_ifft=np.concatenate((rf_specs_ave,rf_specs_ave_conj))
rf_ave=ifft(rf_specs_ifft)
rf_ave_t=np.real(rf_ave)

# %%
fig1,axs1=plt.subplots(nrows=2,ncols=1,figsize=(3,6))
axs1[0].plot(freqs[ii1:ii2],abs(rf_specs_ave[ii1:ii2]),lw=0.6,c='r')
itmax=floor(20/dt)
times=np.linspace(0,(itmax-1)*dt,itmax)
axs1[1].plot(times,rf_ave_t[:itmax],lw=0.6,c='b')
plt.show()

# %%
from scipy.fftpack import ifft,fft
dt=0.1
ttt=np.linspace(0,127*dt,128)
signal=np.sin(2*np.pi*ttt)+np.sin(np.pi*ttt)+np.sin(np.pi/2*ttt)
signal_freq=fft(signal)
signal_recover=ifft(signal_freq)
fig2,axs2=plt.subplots(nrows=2,ncols=1)
axs2[0].plot(abs(signal_freq))
axs2[1].plot(signal)
axs2[1].plot(signal_recover,c='r',lw=0.4)
plt.show()
#%%