from tkinter import *
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from scipy import signal
import math
from tkinter import ttk
import tkinter as tk
import Tools

window=Tk()
window.title("Filtre RII")
window.configure(background = "blue")
window.minsize(width=480,height=400)

def Tracer():
      global x
      freq=fds_val.get()
      x=freq.split(",")
      textbox4.selection_clear()
     
      max=0
      for i in x:
       global amplitude
       global noise 
       global total
       global time
       global total_filt,samplingFrequency
       #print(amplitude)
       if max < float(i):
        max=float(i)
      samplingFrequency= float(fdech_val.get())
      
      samplingInterval =1 / samplingFrequency
      # Begin time period of the signals
      # Time points
      time        = np.arange(beginTime, endTime, samplingInterval)
      amplitude=0

      for i in x:
        amplitude=amplitude + np.sin(2*np.pi*float(i)*time)


      noisepower=float(Pdb_Val.get())
      noise=np.random.normal(0,noisepower,np.size(amplitude))
      total=amplitude+noise

      fch1=filterchoice.get()
      fty1=filtertype.get()
      if fch1=="Filtre Pass Bas":
            f_pass = float(fdarr1_val.get())
            # stop band frequency
            f_stop = float(fdpass1_val.get())
            # pass band ripple
            fs = 0.5
            # pass band freq normalized
            wp = f_pass/(samplingFrequency/2)
            # stop band freq normalized
            ws = f_stop/(samplingFrequency/2)
            # Sampling Time
            Td = 1

            # Conversion to prewrapped analog frequency
            omega_p = (2/Td)*np.tan(wp/2)
            omega_s = (2/Td)*np.tan(ws/2)
            global w,h
            if fty1=="Butterworth":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.buttord(omega_p, omega_s, g_pass, g_stop, analog=True) 
                    if ordre_min.get()==0:
                        N=float(odf_val.get())
                        Wn=wp
                    b, a = signal.butter(N, Wn, 'low', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.filtfilt(b,a,total)
            elif fty1=="Tchebychev type 1":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())                 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb1ord(omega_p, omega_s, g_pass, g_stop, analog=True) 

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.cheby1(N,g_pass, Wn, 'low', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Tchebychev type 2":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())                 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb2ord(omega_p, omega_s, g_pass, g_stop, analog=True)

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.cheby2(N,g_stop, Wn, 'low', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Elliptic":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())                 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.ellipord(omega_p, omega_s, g_pass, g_stop, analog=True)

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.ellip(N,g_pass,g_stop, Wn, 'low', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)
                    

      elif fch1=="Filtre Pass Haut":
            f_pass = float(fdpass1_val.get())
            # stop band frequency
            f_stop = float(fdarr1_val.get())
            # pass band ripple
            fs = 0.5
            # pass band freq normalized
            wp = f_pass/(samplingFrequency/2)
            # stop band freq normalized
            ws = f_stop/(samplingFrequency/2)
            # Sampling Time
            Td = 1
            # Conversion to prewrapped analog frequency
            omega_p = (2/Td)*np.tan(wp/2)
            omega_s = (2/Td)*np.tan(ws/2)
            
            if fty1=="Butterworth":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())     
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.buttord(omega_p, omega_s, g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.butter(N, Wn, 'high', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)
            elif fty1=="Tchebychev type 1":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())                  
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb1ord(omega_p, omega_s, g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.cheby1(N,g_pass, Wn, 'high', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Tchebychev type 2":

                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get()) 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb2ord(omega_p, omega_s, g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.cheby2(N,g_stop, Wn, 'high', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Elliptic":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get()) 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.ellipord(omega_p, omega_s, g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=wp

                    b, a = signal.ellip(N,g_pass,g_stop, Wn, 'high', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)


      elif fch1=="Filtre Pass Bande":
      
            f_stop1 = float(fdarr1_val.get())
            # stop band frequency
            f_pass1 = float(fdpass1_val.get())
            f_pass2 = float(fdpass2_val.get())
            f_stop2 = float(fdarr2_val.get())

            # pass band ripple
            fs = 0.5
            # pass band freq normalized
            wp1 = f_pass1/(samplingFrequency/2)
            wp2 = f_pass2/(samplingFrequency/2)
            # stop band freq normalized
            ws1 = f_stop1/(samplingFrequency/2)
            ws2 = f_stop2/(samplingFrequency/2)
            # Sampling Time
            Td = 1

            # Conversion to prewrapped analog frequency
            omega_p1 = (2/Td)*np.tan(wp1/2)
            omega_p2 = (2/Td)*np.tan(wp2/2)

            omega_s1 = (2/Td)*np.tan(ws1/2)
            omega_s2 = (2/Td)*np.tan(ws2/2)

            if fty1=="Butterworth":
                    # pass band ripple
                    g_pass = float(Adpass_val.get())
                    # stop band attenuation
                    g_stop = float(Adarr_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.buttord([omega_p1,omega_p2], [omega_s1,omega_s2], g_pass, g_stop, analog=True) 

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.butter(N, [wp1,wp2], 'bandpass', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)
            elif fty1=="Tchebychev type 1":
                    # pass band ripple
                    g_pass = float(Adpass_val.get())
                    # stop band attenuation
                    g_stop = float(Adarr_val.get())  
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb1ord([omega_p1,omega_p2], [omega_s1,omega_s2], g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.cheby1(N,g_pass, [wp1,wp2], 'bandpass', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Tchebychev type 2":
                    # pass band ripple
                    g_pass = float(Adpass_val.get())
                    # stop band attenuation
                    g_stop = float(Adarr_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb2ord([omega_p1,omega_p2],[omega_s1,omega_s2], g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.cheby2(N,g_stop, [wp1,wp2], 'bandpass', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Elliptic":
                    # pass band ripple
                    g_pass = float(Adpass_val.get())
                    # stop band attenuation
                    g_stop = float(Adarr_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.ellipord([omega_p1,omega_p2],[omega_s1,omega_s2], g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.ellip(N,g_pass,g_stop, [wp1,wp2], 'bandpass', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)


      elif fch1=="Filtre Coupe Bande":
            
            f_pass1 = float(fdarr1_val.get())
            # stop band frequency
            f_stop1 = float(fdpass1_val.get())
            f_stop2 = float(fdpass2_val.get())
            f_pass2 = float(fdarr2_val.get())

            # pass band ripple
            fs = 0.5
            # pass band freq normalized
            wp1 = f_pass1/(samplingFrequency/2)
            wp2 = f_pass2/(samplingFrequency/2)
            # stop band freq normalized
            ws1 = f_stop1/(samplingFrequency/2)
            ws2 = f_stop2/(samplingFrequency/2)
            # Sampling Time
            Td = 1


            # Conversion to prewrapped analog frequency
            omega_p1 = (2/Td)*np.tan(wp1/2)
            omega_p2 = (2/Td)*np.tan(wp2/2)

            omega_s1 = (2/Td)*np.tan(ws1/2)
            omega_s2 = (2/Td)*np.tan(ws2/2)

            if fty1=="Butterworth":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.buttord([omega_p1,omega_p2], [omega_s1,omega_s2], g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.butter(N, [wp1,wp2], 'bandstop', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)
            elif fty1=="Tchebychev type 1":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())                 
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb1ord([omega_p1,omega_p2], [omega_s1,omega_s2], g_pass, g_stop, analog=True) 

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.cheby1(N,g_pass, [wp1,wp2], 'bandstop', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Tchebychev type 2":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.cheb2ord([omega_p1,omega_p2],[omega_s1,omega_s2], g_pass, g_stop, analog=True)
 
                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.cheby2(N,g_stop, [wp1,wp2], 'bandstop', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)

            elif fty1=="Elliptic":
                    # pass band ripple
                    g_pass = float(Adarr_val.get())
                    # stop band attenuation
                    g_stop = float(Adpass_val.get())
                    # Design of Filter using signal.buttord function
                    N, Wn = signal.ellipord([omega_p1,omega_p2],[omega_s1,omega_s2], g_pass, g_stop, analog=True) 

                    if ordre_min.get()==0:
                         N=float(odf_val.get())
                         Wn=(wp1+wp2)/2

                    b, a = signal.ellip(N,g_pass,g_stop, [wp1,wp2], 'bandstop', False)
                    z, p = signal.bilinear(b, a, fs) # numerator and denominator
                    # w is the freq in z-domain & h is the magnitude in z-domain
                    #global w , h
                    w, h = signal.freqz(z, p, 512) #frequency response 
                    total_filt = signal.lfilter(b,a,total)



      global fourierTransform,frequencies
      fourierTransform = np.fft.fft(amplitude)/len(amplitude)           # Normalize amplitude
      fourierTransform = fourierTransform[range(int(len(amplitude)/2))] # Exclude sampling frequency
      tpCount     = len(amplitude)
      values      = np.arange(int(tpCount/2))
      timePeriod  = tpCount/samplingFrequency
      frequencies = values/timePeriod
      global fourierTransformt,frequenciest

      fourierTransformt = np.fft.fft(total)/len(total)           # Normalize amplitude
      fourierTransformt = fourierTransformt[range(int(len(total)/2))] # Exclude sampling frequency
      tpCountt     = len(total)
      valuest      = np.arange(int(tpCountt/2))
      timePeriodt  = tpCountt/samplingFrequency
      frequenciest = valuest/timePeriodt

      global fourierTransformtf,frequenciestf

      fourierTransformtf = np.fft.fft(total_filt)/len(total_filt)           # Normalize amplitude
      fourierTransformtf = fourierTransformtf[range(int(len(total_filt)/2))] # Exclude sampling frequency
      tpCounttf     = len(total_filt)
      valuestf      = np.arange(int(tpCounttf/2))
      timePeriodtf  = tpCounttf/samplingFrequency
      frequenciestf = valuestf/timePeriodtf

      # Magnitude Reponse
      figure, axis = plt.subplots(4, 1)
      plt.subplots_adjust(hspace=1)
      axis[0].set_title('Signale Sinusoidale avec multiple frequences.')
      axis[0].plot(time, amplitude)
      axis[0].set_xlabel('Temps')
      axis[0].set_ylabel('Amplitude')
      axis[1].set_title('Signale Bruit')
      axis[1].plot(time, total)
      axis[1].set_xlabel('Temps')
      axis[1].set_ylabel('Amplitude')
      axis[2].set_title('Signale Filtree')
      axis[2].plot(time, total_filt)
      axis[2].set_xlabel('Temps')
      axis[2].set_ylabel('Amplitude')
      axis[3].plot(w, 20*np.log10(abs(h)))
      axis[3].set_xscale('log')
      axis[3].set_title('Reponse de la frequence du filtre')
      axis[3].set_xlabel('Frequence [Hz]')
      axis[3].set_ylabel('Amplitude [dB]')
      axis[3].margins(0, 0.1)
      axis[3].grid(which='both', axis='both')
      axis[3].axvline(100, color    ='red')

      figure, axis = plt.subplots(3, 1)
      plt.subplots_adjust(hspace=1)
      axis[0].set_title('Transforme de Fourier du signale')
      axis[0].plot(frequencies, abs(fourierTransform))
      axis[0].set_xlabel('Frequence')
      axis[0].set_ylabel('Amplitude')
      axis[1].set_title('Transforme de Fourier du signale bruit')
      axis[1].plot(frequenciest, abs(fourierTransformt))
      axis[1].set_xlabel('Frequence')
      axis[1].set_ylabel('Amplitude')
      axis[2].set_title('Transforme de Fourier du signale filtree')
      axis[2].plot(frequenciestf, abs(fourierTransformtf))
      axis[2].set_xlabel('Frequence')
      axis[2].set_ylabel('Amplitude')

      plt.show()
        
                       
def on_field_change(index, value, op):
      fch=filterchoice.get()
      if fch=="Filtre Pass Bas":
            fdarr1.configure(text="Frequence de Passage :")
            fdpass1.configure(text="Frequence d'Arret :")
            Adarr.configure(text="Amplitude de Passage :")
            Adpass.configure(text="Amplitude d'Arret :")
            fdpass2.grid_remove()
            fdarr2.grid_remove()
            textbox7.grid_remove()
            textbox8.grid_remove()


      elif fch=="Filtre Pass Haut":
            fdarr1.configure(text="Frequence d'Arret :")
            fdpass1.configure(text="Frequence de Passage :")
            Adarr.configure(text="Amplitude de Passage :")
            Adpass.configure(text="Amplitude d'Arret :")
            fdpass2.grid_remove()
            fdarr2.grid_remove()
            textbox7.grid_remove()
            textbox8.grid_remove()

      elif fch=="Filtre Pass Bande":
            fdpass2.grid(row=12,column=1)
            fdarr2.grid(row=13,column=1)
            fdarr1.configure(text="Frequence d'Arret 1 :")
            fdpass1.configure(text="Frequence de Passage 1 :")
            Adarr.configure(text="Amplitude d'Arret :")
            Adpass.configure(text="Amplitude de Passage :")
          
            fdpass2.configure(text="Frequence de Passage 2 :")
            fdarr2.configure(text="Frequence d'Arret 2 :")

            textbox7.grid(row=12,column=2)
            textbox8.grid(row=13,column=2)

      elif fch=="Filtre Coupe Bande":
            fdpass2.grid(row=12,column=1)
            fdarr2.grid(row=13,column=1)
            fdarr1.configure(text="Frequence de Passage 1 :")
            fdpass1.configure(text="Frequence d'Arret 1 :")
            Adarr.configure(text="Amplitude de Passage :")
            Adpass.configure(text="Amplitude d'Arret :")
          
            fdpass2.configure(text="Frequence d'Arret 2 :")
            fdarr2.configure(text="Frequence de Passage 2 :")

            textbox7.grid(row=12,column=2)
            textbox8.grid(row=13,column=2)


title=Label(text="Filtre RII",fg = "black",bg = "blue",font = "bold ")
title.grid(row=0,column=2)

Pdb=Label(text="La puissance du bruit :",bg = "blue",fg = "yellow",font = "none 9 bold ")
Pdb.grid(row=2,column=1)

Pdb_Val=StringVar()
textbox1 =Entry(textvariable=Pdb_Val)
textbox1.grid(row=2,column=2)

fds=Label(text="Frequence du signale :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fds.grid(row=1,column=1)

fds_val=StringVar()
textbox2=Entry(textvariable=fds_val)
textbox2.grid(row=1,column=2)

tdf=Label(text="Le type du filtre :",bg = "blue",fg = "yellow",font = "none 9 bold ")
tdf.grid(row=3,column=1)

tdf_val = tk.StringVar()
filtertype = ttk.Combobox(window, width = 17, textvariable = tdf_val) 
  # ajouter une combobox des types des filtres
filtertype['values'] = ('Butterworth','Tchebychev type 1','Tchebychev type 2', 'Elliptic')
filtertype.grid(column = 2, row = 3) 
filtertype.current() 

cdf=Label(text="Choix du filtre :",bg = "blue",fg = "yellow",font = "none 9 bold ")
cdf.grid(row=4,column=1)


cdf_val = tk.StringVar() 
filterchoice = ttk.Combobox(window, width = 17, textvariable = cdf_val) 
filterchoice['values'] = ('Filtre Pass Bas','Filtre Pass Haut','Filtre Pass Bande', 'Filtre Coupe Bande')                 
filterchoice.grid(column = 2, row = 4) 
filterchoice.current() 
cdf_val.trace('w',on_field_change)

odf=Label(text="L'ordre du filtre :",bg = "blue",fg = "yellow",font = "none 9 bold ")
odf.grid(row=8,column=1)


odf_val=StringVar()
textbox3=Entry(textvariable=odf_val)
textbox3.grid(row=8,column=2)

ordre_min = IntVar()
check=Checkbutton(window, text="L'ordre Minimale",bg = "lightgreen", variable=ordre_min)
check.grid(row=8,column=3, sticky=W)

common_imag = PhotoImage(width = 1, height = 1)
button=Button(text="Ploter la Solution",image = common_imag,compound = "c", bg = "lightgreen",command=Tracer)
button.grid(row=20,column=2)

fdech=Label(text="Frequence d'Echantillonage :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fdech.grid(row=9,column=1)

fdech_val=StringVar()
textbox4=Entry(textvariable=fdech_val)
textbox4.grid(row=9,column=2)

fdarr1=Label(text="Frequence d'Arret 1 :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fdarr1.grid(row=10,column=1)

fdarr1_val=StringVar()
textbox5=Entry(textvariable=fdarr1_val)
textbox5.grid(row=10,column=2)

fdpass1=Label(text="Frequence de Passage 1 :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fdpass1.grid(row=11,column=1)

fdpass1_val=StringVar()
textbox6=Entry(textvariable=fdpass1_val)
textbox6.grid(row=11,column=2)

fdpass2=Label(text="Frequence de Passage 2 :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fdpass2.grid(row=12,column=1)

fdpass2_val=StringVar()
textbox7=Entry(textvariable=fdpass2_val)
textbox7.grid(row=12,column=2)

fdarr2=Label(text="Frequence d'Arret 2 :",bg = "blue",fg = "yellow",font = "none 9 bold ")
fdarr2.grid(row=13,column=1)

fdarr2_val=StringVar()
textbox8=Entry(textvariable=fdarr2_val)
textbox8.grid(row=13,column=2)

Adarr=Label(text="Amplitude d'Arret :",bg = "blue",fg = "yellow",font = "none 9 bold ")
Adarr.grid(row=14,column=1)

Adarr_val=StringVar()
textbox9=Entry(textvariable=Adarr_val)
textbox9.grid(row=14,column=2)

Adpass=Label(text="Amplitude de Passage :",bg = "blue",fg = "yellow",font = "none 9 bold ")
Adpass.grid(row=15,column=1)

Adpass_val=StringVar()
textbox10=Entry(textvariable=Adpass_val)
textbox10.grid(row=15,column=2)


amplitude=0
beginTime           = 0
      # End time period of the signals
endTime             = 10


window.mainloop()
