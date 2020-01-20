#!/usr/bin/env python
#
# Paul's Extreme Sound Stretch (Paulstretch) - Python version
# using a new method
#
# by Nasca Octavian PAUL, Targu Mures, Romania
# http://www.paulnasca.com/
#
# http://hypermammut.sourceforge.net/paulstretch/
#
# this file is released under Public Domain
#


import sys
from numpy import *
import scipy.io.wavfile
import wave
from optparse import OptionParser

plot_onsets=False
if plot_onsets:
    import matplotlib.pyplot as plt


def load_wav(filename):
    try:
        wavedata=scipy.io.wavfile.read(filename)
        samplerate=int(wavedata[0])
        smp=wavedata[1]*(1.0/32768.0)
        smp=smp.transpose()
        if len(smp.shape)==1: #convert to stereo
            smp=tile(smp,(2,1))
        return (samplerate,smp)
    except:
        print ("Error loading wav: "+filename)
        return None



def optimize_windowsize(n):
    orig_n=n
    while True:
        n=orig_n
        while (n%2)==0:
            n/=2
        while (n%3)==0:
            n/=3
        while (n%5)==0:
            n/=5

        if n<2:
            break
        orig_n+=1
    return orig_n

def paulstretch(samplerate,smp,stretch,windowsize_seconds,onset_level,outfilename):

    if plot_onsets:
        onsets=[]

    nchannels=smp.shape[0]

    outfile=wave.open(outfilename,"wb")
    outfile.setsampwidth(2)
    outfile.setframerate(samplerate)
    outfile.setnchannels(nchannels)

    #make sure that windowsize is even and larger than 16
    windowsize=int(windowsize_seconds*samplerate)
    if windowsize<16:
        windowsize=16
    windowsize=optimize_windowsize(windowsize)
    windowsize=int(windowsize/2)*2
    half_windowsize=int(windowsize/2)

    #correct the end of the smp
    nsamples=smp.shape[1]
    end_size=int(samplerate*0.05)
    if end_size<16:
        end_size=16

    smp[:,nsamples-end_size:nsamples]*=linspace(1,0,end_size)

    
    #compute the displacement inside the input file
    start_pos=0.0
    displace_pos=windowsize*0.5

    #create Hann window
    window=0.5-cos(arange(windowsize,dtype='float')*2.0*pi/(windowsize-1))*0.5

    old_windowed_buf=zeros((2,windowsize))
    hinv_sqrt2=(1+sqrt(0.5))*0.5
    hinv_buf=2.0*(hinv_sqrt2-(1.0-hinv_sqrt2)*cos(arange(half_windowsize,dtype='float')*2.0*pi/half_windowsize))/hinv_sqrt2

    freqs=zeros((2,half_windowsize+1))
    old_freqs=freqs

    num_bins_scaled_freq=32
    freqs_scaled=zeros(num_bins_scaled_freq)
    old_freqs_scaled=freqs_scaled

    displace_tick=0.0
    displace_tick_increase=1.0/stretch
    if displace_tick_increase>1.0:
        displace_tick_increase=1.0
    extra_onset_time_credit=0.0
    get_next_buf=True
    while True:
        if get_next_buf:
            old_freqs=freqs
            old_freqs_scaled=freqs_scaled

            #get the windowed buffer
            istart_pos=int(floor(start_pos))
            buf=smp[:,istart_pos:istart_pos+windowsize]
            if buf.shape[1]<windowsize:
                buf=append(buf,zeros((2,windowsize-buf.shape[1])),1)
            buf=buf*window
    
            #get the amplitudes of the frequency components and discard the phases
            freqs=abs(fft.rfft(buf))

            #scale down the spectrum to detect onsets
            freqs_len=freqs.shape[1]
            if num_bins_scaled_freq<freqs_len:
                freqs_len_div=freqs_len//num_bins_scaled_freq
                new_freqs_len=freqs_len_div*num_bins_scaled_freq
                freqs_scaled=mean(mean(freqs,0)[:new_freqs_len].reshape([num_bins_scaled_freq,freqs_len_div]),1)
            else:
                freqs_scaled=zeros(num_bins_scaled_freq)


            #process onsets
            m=2.0*mean(freqs_scaled-old_freqs_scaled)/(mean(abs(old_freqs_scaled))+1e-3)
            if m<0.0:
                m=0.0
            if m>1.0:
                m=1.0
            if plot_onsets:
                onsets.append(m)
            if m>onset_level:
                displace_tick=1.0
                extra_onset_time_credit+=1.0

        cfreqs=(freqs*displace_tick)+(old_freqs*(1.0-displace_tick))

        #randomize the phases by multiplication with a random complex number with modulus=1
        ph=random.uniform(0,2*pi,(nchannels,cfreqs.shape[1]))*1j
        cfreqs=cfreqs*exp(ph)

        #do the inverse FFT 
        buf=fft.irfft(cfreqs)

        #window again the output buffer
        buf*=window

        #overlap-add the output
        output=buf[:,0:half_windowsize]+old_windowed_buf[:,half_windowsize:windowsize]
        old_windowed_buf=buf

        #remove the resulted amplitude modulation
        output*=hinv_buf
        
        #clamp the values to -1..1 
        output[output>1.0]=1.0
        output[output<-1.0]=-1.0

        #write the output to wav file
        outfile.writeframes(int16(output.ravel('F')*32767.0).tostring())

        if get_next_buf:
            start_pos+=displace_pos

        get_next_buf=False

        if start_pos>=nsamples:
            print ("100 %")
            break
        sys.stdout.write ("%d %% \r" % int(100.0*start_pos/nsamples))
        sys.stdout.flush()

        
        if extra_onset_time_credit<=0.0:
            displace_tick+=displace_tick_increase
        else:
            credit_get=0.5*displace_tick_increase #this must be less than displace_tick_increase
            extra_onset_time_credit-=credit_get
            if extra_onset_time_credit<0:
                extra_onset_time_credit=0
            displace_tick+=displace_tick_increase-credit_get

        if displace_tick>=1.0:
            displace_tick=displace_tick % 1.0
            get_next_buf=True

    outfile.close()
    
    if plot_onsets:
        plt.plot(onsets)
        plt.show()
    

########################################
print ("Paul's Extreme Sound Stretch (Paulstretch) - Python version 20141220")
print ("new method: using onsets information")
print ("by Nasca Octavian PAUL, Targu Mures, Romania\n")
parser = OptionParser(usage="usage: %prog [options] input_wav output_wav")
parser.add_option("-s", "--stretch", dest="stretch",help="stretch amount (1.0 = no stretch)",type="float",default=8.0)
parser.add_option("-w", "--window_size", dest="window_size",help="window size (seconds)",type="float",default=0.25)
parser.add_option("-t", "--onset", dest="onset",help="onset sensitivity (0.0=max,1.0=min)",type="float",default=10.0)
(options, args) = parser.parse_args()


if (len(args)<2) or (options.stretch<=0.0) or (options.window_size<=0.001):
    print ("Error in command line parameters. Run this program with --help for help.")
    sys.exit(1)

print ("stretch amount = %g" % options.stretch)
print ("window size = %g seconds" % options.window_size)
print ("onset sensitivity = %g" % options.onset)
(samplerate,smp)=load_wav(args[0])

paulstretch(samplerate,smp,options.stretch,options.window_size,options.onset,args[1])



