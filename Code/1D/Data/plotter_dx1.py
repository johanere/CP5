
import numpy as np
import matplotlib.pyplot as plt
import string as str
import matplotlib.pyplot as plt

def tex(x):
    exp = int(np.floor(np.log10(abs(x))))
    number=x/10**exp
    return   "$%.2g \cdot 10^{%s}$"%(number,exp)


outdata=[[],[],[],[]]
outdata[0]=["results_FW_0.05T_0.1dx","results_BW_0.05T_0.1dx","results_CN_0.05T_0.1dx"]
outdata[1]=["results_FW_0.2T_0.1dx","results_BW_0.2T_0.1dx","results_CN_0.2T_0.1dx"]
outdata[2]=["results_FW_0.05T_0.01dx","results_BW_0.05T_0.01dx","results_CN_0.05T_0.01dx"]
outdata[3]=["results_FW_0.2T_0.01dx","results_BW_0.2T_0.01dx","results_CN_0.2T_0.01dx"]

figurtitle=["../../Results/Comparison_1.png", "../../Results/Comparison_2.png"]
tab=np.zeros((3,2))

average_error=np.zeros((3,4))
max_error=np.zeros((3,4))



#open for reading
for step in range(0,2):
    f, axarr = plt.subplots(1, 2, sharex='all',figsize=(9, 4))
    for l in range(0,2):
        for j in range(len(outdata[l+step*2])):
            f= open(outdata[l+step*2][j])
            if f.mode == 'r':
                lines = f.readlines()
            N=len(lines)
    
            method=lines[0].strip()
            gridsize=lines[1].split()[0]
            n=lines[2].split()[0]
            T=float(lines[3].split()[0])
            dx=float(lines[4].split()[0])
            dt=lines[5].split()[0]
            
            x=np.linspace(0,1,N-7)
            
            solution=np.zeros(N-7)
            exact=np.zeros(N-7)
            relerr=np.zeros(N-7)
            
            for i in range(0,N-7):
                hold=lines[i+7].split()
                solution[i]=float(hold[0])
                exact[i]=float(hold[1])
                relerr[i]=float(hold[2])
            
            max_error[j,l+step*2]=max(relerr)
            average_error[j,l+step*2]=sum(relerr)/len(relerr)

    
            axarr[l].plot(x,solution,label=method)
            

        
        axarr[l].plot(x,exact,label="Analytical")
        axarr[l].set_title("$t_f=%2g$"%(T),fontsize=10)
        
    for i in range (2):
        axarr[i].axis('equal')
        axarr[i].axis('equal')
        axarr[i].grid(True, linestyle='-.')
        axarr[i].legend()
        axarr[i].set_xlabel('$x$')
        axarr[i].set_ylabel('$temperature$')
    plt.savefig("%s"%figurtitle[step], bbox_inches='tight')

#plt.show()

outfile=open("../../Results/tab_maxerror.txt","w") 
for i in range(0,3):
    for j in range(0,4):
        outfile.write(tex(max_error[i,j]))
        if j==3:
            outfile.write("\\\ \n")
        else:
            outfile.write(" & ")  
outfile.close()

outfile=open("../../Results/tab_average.txt","w") 
for i in range(0,3):
    for j in range(0,4):
        outfile.write(tex(average_error[i,j]))
        if j==3:
            outfile.write("\\\ \n")
        else:
            outfile.write(" & ")  
outfile.close()

