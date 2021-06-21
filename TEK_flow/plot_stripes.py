import matplotlib.pyplot as plt
import tikzplotlib
import numpy as np



s1 = 'maroon'
s2 = 'firebrick'
s3 = 'salmon'
s4 = 'sandybrown'
s5 = 'orange'
s6 = 'gold'
s7 = 'yellow'
s8 = 'greenyellow'
s9 = 'seagreen'
s10 = 'skyblue'
s11 = 'cornflowerblue'
s12 = 'slateblue'
s13 = 'darkblue'
s14 = 'blue' 

def PLOT(label,Cmin,Cmax,color,filename,t0=None):
    time = np.loadtxt(filename,comments="#",usecols=0)[Cmin:Cmax]
    data = np.loadtxt(filename,comments="#",usecols=1)[Cmin:Cmax]
    jerr = np.loadtxt(filename,comments="#",usecols=2)[Cmin:Cmax]

    if not t0==None: time = time/t0

    plt.fill_between(time,data+jerr,data-jerr,color=color,label=label,alpha=0.7)  # se vuoi shade


plt.figure(figsize=(3.5,7))
plt.ylabel(r'$\frac{<t^2E(t)>}{N}$')
plt.xlabel(r'$\frac{t}{a^2}$')
plt.grid(zorder=0, linewidth=0.5, color='gainsboro')

# ==================================================================
# ------------------------ nfhN289b036 -----------------------------
'''
title = 'nfhn289b036'
plt.title(r'$N_f=\frac{1}{2}$, $N=289$, $b=0.36$')
PLOT(r'$\kappa_a=0.100$',0,320,s1,"nfh/n289b0360k5hf100.dat")
PLOT(r'$\kappa_a=0.110$',0,320,s2,"nfh/n289b0360k5hf110.dat")
PLOT(r'$\kappa_a=0.110$',0,320,s3,"nfh/n289b0360k5hf120.dat")
PLOT(r'$\kappa_a=0.130$',0,320,s4,"nfh/n289b0360k5hf130.dat")
PLOT(r'$\kappa_a=0.140$',0,320,s5,"nfh/n289b0360k5hf140.dat")
PLOT(r'$\kappa_a=0.150$',0,320,s6,"nfh/n289b0360k5hf150.dat")
PLOT(r'$\kappa_a=0.155$',0,320,s7,"nfh/n289b0360k5hf155.dat")
PLOT(r'$\kappa_a=0.160$',0,320,s8,"nfh/n289b0360k5hf160.dat")
PLOT(r'$\kappa_a=0.165$',0,320,s9,"nfh/n289b0360k5hf165.dat")
PLOT(r'$\kappa_a=0.175$',0,320,s10,"nfh/n289b0360k5hf175.dat")
'''
# --------------------------- nfhN289b035 ---------------------------
'''
title = 'nfhn289b035_norm'
plt.title(r'$N=289$, $b=0.35$')#'$N_f=\frac{1}{2}$, $N=289$, $b=0.35$')
PLOT(r"$\kappa=0.100$",0,320,s1,"nfh/n289b0350k5hf100.dat")
PLOT(r"$\kappa=0.110$",0,320,s2,"nfh/n289b0350k5hf110.dat")
PLOT(r"$\kappa=0.110$",0,320,s3,"nfh/n289b0350k5hf120.dat")
PLOT(r"$\kappa=0.130$",0,320,s4,"nfh/n289b0350k5hf130.dat")
PLOT(r"$\kappa=0.140$",0,320,s5,"nfh/n289b0350k5hf140.dat")
PLOT(r"$\kappa=0.150$",0,320,s6,"nfh/n289b0350k5hf150.dat")
PLOT(r"$\kappa=0.155$",0,320,s7,"nfh/n289b0350k5hf155.dat")
PLOT(r"$\kappa=0.160$",0,320,s8,"nfh/n289b0350k5hf160.dat")
PLOT(r"$\kappa=0.165$",0,320,s9,"nfh/n289b0350k5hf170.dat")
PLOT(r"$\kappa=0.175$",0,320,s10,"nfh/n289b0350k5hf175.dat")
PLOT(r"$\kappa_{adj}=0.1775$",0,640,s1,"nfh/n289b0350k5hf1775_norm.dat")
PLOT(r"$\kappa_{adj}=0.180$", 0,640,s2,"nfh/n289b0350k5hf1800_norm.dat")
PLOT(r"$\kappa_{adj}=0.1825$",0,640,s3,"nfh/n289b0350k5hf1825_norm.dat")
PLOT(r"$\kappa_{adj}=0.185$", 0,640,s4,"nfh/n289b0350k5hf1850_norm.dat")
PLOT(r"$\kappa_{adj}=0.1875$",0,640,s5,"nfh/n289b0350k5hf1875_norm.dat")
'''
# ------------------------ nfhN169b035 -----------------------------
'''
title = 'nfhn169b035_norm'
plt.title(r'$N=169$, $b=0.35$')#'$N_f=\frac{1}{2}$, $N=289$, $b=0.35$')
PLOT(r'$\kappa_{adj}=0.1775$',0,640,s1,"nfh/n169b0350k5hf1775_norm.dat")
PLOT(r'$\kappa_{adj}=0.1800$',0,640,s2,"nfh/n169b0350k5hf1800_norm.dat")
PLOT(r'$\kappa_{adj}=0.1825$',0,640,s3,"nfh/n169b0350k5hf1825_norm.dat")
PLOT(r'$\kappa_{adj}=0.1850$',0,640,s4,"nfh/n169b0350k5hf1850_norm.dat")
'''
'''
title = 'nfhn169b034_norm'
plt.title(r'$N=169$, $b=0.34$')#'$N_f=\frac{1}{2}$, $N=289$, $b=0.35$')
PLOT(r'$\kappa_{adj}=0.1875$',0,640,s1,"nfh/n169b0340k5hf1875_norm.dat")
PLOT(r'$\kappa_{adj}=0.1890$',0,640,s2,"nfh/n169b0340k5hf1890_norm.dat")
PLOT(r'$\kappa_{adj}=0.1910$',0,640,s3,"nfh/n169b0340k5hf1910_norm.dat")
'''
'''
title = 'nfhn169b0347_norm'
plt.title(r'$N=169$, $b=0.347$')#'$N_f=\frac{1}{2}$, $N=289$, $b=0.35$')
PLOT(r'$\kappa_{adj}=0.1850$',0,640,s1,"nfh/n169b0347k5hf1850_norm.dat")
PLOT(r'$\kappa_{adj}=0.1875$',0,640,s2,"nfh/n169b0347k5hf1875_norm.dat")
PLOT(r'$\kappa_{adj}=0.1890$',0,640,s3,"nfh/n169b0347k5hf1890_norm.dat")
'''
# ------------------------ nfhN361b035 -----------------------------
title = 'nfhn361b035_norm'
plt.title(r'$N=361$, $b=0.35$')#'$N_f=\frac{1}{2}$, $N=289$, $b=0.35$')
PLOT(r'$\kappa_{adj}=0.1825$',0,640,s1,"DATA/nfh/n361b0350k7hf1825_norm.dat")
PLOT(r'$\kappa_{adj}=0.1875$',0,640,s4,"DATA/nfh/n361b0350k7hf1875_norm.dat")
'''
'''
# ======================================================



# ================== nf1b036 ============================
'''
plt.title(r'$N_f=1$, N=121, b=0.36')
PLOT(r'$\kappa=0.120$',0,320,s3,"nf1/n121b0360k3of120.dat")
PLOT(r'$\kappa=0.130$',0,320,s4,"nf1/n121b0360k3of130.dat")
PLOT(r'$\kappa=0.140$',0,320,s5,"nf1/n121b0360k3of140.dat")
PLOT(r'$\kappa=0.150$',0,320,s6,"nf1/n121b0360k3of150.dat")
PLOT(r'$\kappa=0.160$',0,320,s8,"nf1/n121b0360k3of160.dat")
PLOT(r'$\kappa=0.165$',0,320,s9,"nf1/n121b0360k3of165.dat")
PLOT(r'$\kappa=0.170$',0,320,s10,"nf1/n121b0360k3of170.dat")
'''
# ------------------------ nf1b035 ---------------
'''
title = 'nf1n121b035'
plt.title(r'$N_f=1$, $N=121$, $b=0.35$')
#PLOT(r'$\kappa=0.120$',0,320,s1,"nf1/n121b0350k3of120.dat")
PLOT(r'$\kappa=0.130$',0,320,s2,"nf1/n121b0350k3of130.dat")
PLOT(r'$\kappa=0.140$',0,320,s3,"nf1/n121b0350k3of140.dat")
PLOT(r'$\kappa=0.150$',0,320,s4,"nf1/n121b0350k3of150.dat")
PLOT(r'$\kappa=0.155$',0,320,s5,"nf1/n121b0350k3of155.dat")
'''
'''
title = 'nf1n289b035'
plt.title(r'$N_f=1$, $N=289$, $b=0.35$')
PLOT(r'$\kappa=0.130$',0,320,s2,"nf1/n289b0350k5of130.dat",2.32226454571327)
PLOT(r'$\kappa=0.140$',0,320,s3,"nf1/n289b0350k5of140.dat",2.67251796593814)
PLOT(r'$\kappa=0.150$',0,320,s4,"nf1/n289b0350k5of150.dat",3.42004836956662)
PLOT(r'$\kappa=0.160$',0,320,s6,"nf1/n289b0350k5of160.dat",5.31423573765640)
PLOT(r'$\kappa=0.170$',0,320,s7,"nf1/n289b0350k5of170.dat",10.2614910639642)
'''
'''
title = 'nf1n289b036'
plt.title(r'$N_f=1$, $N=289$, $b=0.36$')
PLOT(r'$\kappa=0.120$',0,320,s1,"nf1/n289b0360k5of120.dat")
PLOT(r'$\kappa=0.130$',0,320,s2,"nf1/n289b0360k5of130.dat")
PLOT(r'$\kappa=0.140$',0,320,s3,"nf1/n289b0360k5of140.dat")
PLOT(r'$\kappa=0.150$',0,320,s4,"nf1/n289b0360k5of150.dat")
PLOT(r'$\kappa=0.170$',0,320,s5,"nf1/n289b0360k5of155.dat")
PLOT(r'$\kappa=0.160$',0,320,s6,"nf1/n289b0360k5of160.dat")
PLOT(r'$\kappa=0.160$',0,320,s7,"nf1/n289b0360k5of165.dat")
'''
# ======================================================


# ================== nf2b035 ============================
'''
title = 'nf2n289b035'
plt.title(r'$N_f=2$, $N=289$, $b=0.35$')
PLOT(r'$\kappa=0.10$',0,320,s1,"nf2/n289b0350k5kp100.dat",2.0128823399456466)
PLOT(r'$\kappa=0.11$',0,320,s2,"nf2/n289b0350k5kp110.dat",2.2183042021525280)
PLOT(r'$\kappa=0.12$',0,320,s3,"nf2/n289b0350k5kp120.dat",2.5666218941968410)
PLOT(r'$\kappa=0.13$',0,320,s4,"nf2/n289b0350k5kp130.dat",3.2022665336523928)
PLOT(r'$\kappa=0.14$',0,320,s5,"nf2/n289b0350k5kp140.dat",4.5616691899598700)
PLOT(r'$\kappa=0.15$',0,320,s6,"nf2/n289b0350k5kp150.dat",9.3375281530120900)
#PLOT(r'$\kappa=0.16$',0,320,s7,"nf2/n289b0350k5kp160.dat",)
#PLOT(r'$\kappa=0.17$',0,320,s8,"nf2/n289b0350k5kp170.dat",)
'''
'''
title = 'nf2n289b036'
plt.title(r'$N_f=2$, $N=289$, $b=0.36$')
PLOT(r'$\kappa=0.10$',0,320,s1,"nf2/n289b0360k5kp100.dat")
PLOT(r'$\kappa=0.11$',0,320,s2,"nf2/n289b0360k5kp110.dat")
PLOT(r'$\kappa=0.12$',0,320,s3,"nf2/n289b0360k5kp120.dat")
PLOT(r'$\kappa=0.13$',0,320,s4,"nf2/n289b0360k5kp130.dat")
PLOT(r'$\kappa=0.14$',0,320,s5,"nf2/n289b0360k5kp140.dat")
PLOT(r'$\kappa=0.15$',0,320,s6,"nf2/n289b0360k5kp150.dat")
PLOT(r'$\kappa=0.16$',0,320,s7,"nf2/n289b0360k5kp160.dat")
PLOT(r'$\kappa=0.17$',0,320,s8,"nf2/n289b0360k5kp170.dat")
'''
# ------------------------ nf2b035 ---------------
'''
title = 'nf2n121b035'
plt.title(r'$N_f=2$, $N=121$, $b=0.35$')
PLOT(r'$\kappa=0.120$',0,320,s3,"nf2/n121b0350k3kp120.dat")
PLOT(r'$\kappa=0.130$',0,320,s4,"nf2/n121b0350k3kp130.dat")
PLOT(r'$\kappa=0.140$',0,320,s5,"nf2/n121b0350k3kp140.dat")
PLOT(r'$\kappa=0.150$',0,320,s6,"nf2/n121b0350k3kp150.dat")
PLOT(r'$\kappa=0.160$',0,320,s7,"nf2/n121b0350k3kp160.dat")
PLOT(r'$\kappa=0.165$',0,320,s8,"nf2/n121b0350k3kp165.dat")
PLOT(r'$\kappa=0.170$',0,320,s9,"nf2/n121b0350k3kp170.dat")
PLOT(r'$\kappa=0.175$',0,320,s10,"nf2/n121b0350k3kp175.dat")
'''
'''
title = 'nf2n121b036'
plt.title(r'$N_f=2$, $N=121$, $b=0.36$')
PLOT(r'$\kappa=0.120$',0,320,s3,"nf2/n121b0360k3kp120.dat")
PLOT(r'$\kappa=0.130$',0,320,s4,"nf2/n121b0360k3kp130.dat")
PLOT(r'$\kappa=0.140$',0,320,s5,"nf2/n121b0360k3kp140.dat")
PLOT(r'$\kappa=0.150$',0,320,s6,"nf2/n121b0360k3kp150.dat")
PLOT(r'$\kappa=0.160$',0,320,s7,"nf2/n121b0360k3kp160.dat")
PLOT(r'$\kappa=0.165$',0,320,s8,"nf2/n121b0360k3kp165.dat")
PLOT(r'$\kappa=0.170$',0,320,s9,"nf2/n121b0360k3kp170.dat")
PLOT(r'$\kappa=0.175$',0,320,s10,"nf2/n121b0360k3kp175.dat")
'''
# =======================================================

plt.show()
