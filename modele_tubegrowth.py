
# Parametres du modele
import random
import math
import pandas as pd

bed="Parameter.txt" 
bed = open(str(bed), "r")
bed=bed.read()
bed = bed.split("#")

N=bed[1]
N = N.split("\n")
N=int(N[1])
individus=range(1,N+1)

nb_generation=bed[2]
nb_generation = nb_generation.split("\n")
nb_generation=int(nb_generation[1])

replica=bed[3]
replica = replica.split("\n")
replica=int(replica[1])

mutation_rate=["y"]
mutation_rate2=bed[4]
mutation_rate2=mutation_rate2.split("\n")
mutation_rate2=int(mutation_rate2[1])
while len(mutation_rate)<mutation_rate2:
    mutation_rate.append("n")


## genes

###pollen tube growth V
speed=range(0,11)
SDspeed=3

### style length
SL=range(0,11)
SDSL=3
OVUL=10


## tested
autofec=bed[5]
autofec = autofec.split("\n")
autofec=float(autofec[1])

ratio=bed[6]
ratio=ratio.split("\n")
ratio=int(ratio[1])

cost=bed[7]
cost=cost.split("\n")
cost=float(cost[1])

donnor=bed[8]
donnor=donnor.split("\n")
donnor=int(donnor[1])

intra_selection=bed[9]
intra_selection=intra_selection.split("\n")
intra_selection=float(intra_selection[1])

coef_survie_stigmat=bed[10]
coef_survie_stigmat=coef_survie_stigmat.split("\n")
coef_survie_stigmat=float(coef_survie_stigmat[1])

coef_survie_pollen=bed[11]
coef_survie_pollen=coef_survie_pollen.split("\n")
coef_survie_pollen=float(coef_survie_pollen[1])


#sorties modele coevulution non lié


fichier = open("tube growth_positive_NO_link_N_"+str(N)+"_selfing_"+str(autofec)+"_ratio_"+str(ratio)+"_female_cost_"+str(cost)+"_donnor_"+str(donnor)+"_femelle_"+str(intra_selection)+"_proxy_femelle_"+str(coef_survie_stigmat)+"_proxy_male_"+str(coef_survie_pollen)+".csv", "a")
fichier.write("simulation;generation;mean_speed_pt;mean_size_style;mean_speed_neutral;mean_size_neutral;nb_graine")
fichier2 = open("tube growth_positive_NO_link_N_"+str(N)+"_selfing_"+str(autofec)+"_ratio_"+str(ratio)+"_female_cost_"+str(cost)+"_donnor_"+str(donnor)+"_femelle_"+str(intra_selection)+"_proxy_femelle_"+str(coef_survie_stigmat)+"_proxy_male_"+str(coef_survie_pollen)+"HAPLOTYPE.csv", "a")
fichier2.write("simulation;generation;ind;")

l=1
while l <= replica :
    ind=0
    pop_initiale_1=[]
    #creation des genotypes de chaque individu de la pop 
    while ind<=len(individus)-1:
        ptg1=random.choice(speed)
        x1=ptg1
        if x1<0:x1=0
        ptg1=random.choice(speed)
        x2=ptg1
        if x2<0:x2=0
        x5=random.choice(speed)
        if x5<0:x5=0
        x6=random.choice(speed)
        if x6<0:x6=0
        ptg1=random.choice(SL)
        x3=ptg1
        if x3<0:x3=0
        ptg1=random.choice(SL)
        x4=ptg1
        if x4<0:x4=0
        x7=random.choice(SL)
        x8=random.choice(SL)
        if x7<0:x7=0
        if x8<0:x8=0
        #taille de pistil
        size=x3+x4
        #nombre d'ovule
        OVUL2=round(OVUL-cost*size,0)
        if OVUL2<=0:OVUL2=0
        # nomdre d'ovule survivants
        survi=1
        # nombre de pollen acceuilli
        partner=ratio*OVUL
        genotype_S_1=[x1,x2,x3,x4,x5,x6,x7,x8,size,OVUL2,survi,partner] #OK
        pop_initiale_1.append(genotype_S_1)
        ind=ind+1
    #sortie modele
    x="\n"+str(l)+";0"
    i=0
    SPEED=0
    STYLE=0
    spe_neu=0
    sty_neu=0
    while i <= len(pop_initiale_1)-1:
        ind=pop_initiale_1[i]
        fichier2.write(str(x)+";"+str(ind))
        SPEED=SPEED+ind[0]+ind[1]
        STYLE=STYLE+ind[2]+ind[3]
        spe_neu=spe_neu+ind[4]+ind[5]
        sty_neu=sty_neu+ind[6]+ind[7]
        i=i+1
    SPEED=SPEED/(2*N)
    STYLE=STYLE/(2*N)
    spe_neu=spe_neu/(2*N)
    sty_neu=sty_neu/(2*N)
    x=str(x)+";"+str(SPEED)+";"+str(STYLE)+";"+str(spe_neu)+";"+str(sty_neu)+";NA"#ok
    fichier.write(str(x))
    # faire les croisements
    generation_ini=1
    suivi1=pop_initiale_1
    ## formation des descendants
    while generation_ini<=nb_generation:
        creation1=[]
        table1=suivi1[:]
        random.shuffle(table1)
        ind=0
        graine=[]
        while ind<N and ind<len(table1)-1:
            mere=table1[ind]
            #formation des ovules
            ovule=[]
            if mere[9]>0:
                while len(ovule)<mere[9]:
                    x=[random.choice([mere[0],mere[1]])]
                    X=random.choice([mere[2],mere[3]])
                    x.append(X)
                    X=random.choice([mere[4],mere[5]])
                    x.append(X)
                    X=random.choice([mere[6],mere[7]])
                    x.append(X)
                    ovule.append(x)
                # tirage des pollens donnor  sur stigmate
                table2=[]
                i=0
                while i<=len(table1)-1:
                    if i!=ind:
                        table2.append(table1[i])
                    i=i+1
                random.shuffle(table2)
                if autofec==0: table2=table2[0:donnor]
                elif autofec!=0: table2=table2[0:donnor-1]
                                   
                pollen=[]

                while len(pollen) < (autofec*ratio*mere[9]):
                    x=[random.choice([mere[0],mere[1]])]
                    X=random.choice([mere[2],mere[3]])
                    x.append(X)
                    X=random.choice([mere[4],mere[5]])
                    x.append(X)
                    X=random.choice([mere[6],mere[7]])
                    x.append(X)
                    pollen.append(x)
                while len(pollen)<mere[11]:
                    xb=random.choice(table2)
                    x=[random.choice([xb[0],xb[1]])]
                    X=random.choice([xb[2],xb[3]])
                    x.append(X)
                    X=random.choice([xb[4],xb[5]])
                    x.append(X)
                    X=random.choice([xb[6],xb[7]])
                    x.append(X)
                    pollen.append(x)
                
                pollen = pd.DataFrame(pollen, index = range(0,len(pollen)), columns = ['V', 'L', 'NtV', 'NTL'])
                # nb de pollen pouvant féconder les ovules (choix intersex)
                if float(mere[8])!=0: intra=int((ratio*OVUL)/((float(mere[8]))*intra_selection))
                if float(mere[8])==0: intra=len(pollen)
                
                # taille des pollens
                pollen=pollen.sort_values(by=["V"], ascending = False)
                if intra<len(pollen.V):
                    seq=pollen.iloc[intra-1]
                    seq=seq[0]
                    ptg2=pollen[(pollen.V>=float(seq))]
                elif intra>=len(pollen.V):
                    ptg2=pollen
                ptg2=ptg2.values.tolist()
                random.shuffle(ptg2)
                
                # fécondation
                ptg2=ptg2[0:int(mere[9])]
                ptg=ptg2
                i=0
                random.shuffle(ptg)
                x=ptg
                random.shuffle(ovule)
                x2=ovule
                while i<=int(mere[9])-1 and i<=len(ptg)-1:
                    G1=x[i]
                    G2=x2[i]
                    G=[G1[0],G2[0],G1[1],G2[1],G1[2],G2[2],G1[3],G2[3]]
                    graine.append(G)
                    i=i+1

            ind=ind+1
        
        # proxy survie
        random.shuffle(graine)
        if len(graine)<=N :graine2=graine[:]
        
        elif  len(graine)>N and coef_survie_stigmat!=0 or len(graine)>N and coef_survie_pollen!=0:
            graine = pd.DataFrame(graine, index = range(0,len(graine)), columns = ['V1','V2','L1','L2','NtV1','NtV2','NTL1','NTL2'])
            graine["G1"]=1+(coef_survie_stigmat*((graine.L1+graine.L2)))+(coef_survie_pollen*((graine.V1+graine.V2)))
            graine=graine.sort_values(by=["G1"], ascending = False)
            seq=graine.iloc[N]
            seq=seq[8]
            graine2=graine[(graine.G1>=float(seq))]
            graine2=graine2.values.tolist()
            
            

        elif len(graine)>N and coef_survie_stigmat==0 and coef_survie_pollen==0:
            random.shuffle(graine)
            graine2=graine[0:N]
            



                
            
        x=graine2
        i=0
        while i<=len(x)-1 and len(creation1)<N:
            G1=x[i]
            ptg1=G1[0]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[0])-SDspeed,int(G1[0])+SDspeed)
                ptg1=random.choice(ptg)
            x1=ptg1
            if x1<0:x1=0
            ptg1=G1[1]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[1])-SDspeed,int(G1[1])+SDspeed)
                ptg1=random.choice(ptg)
            x2=ptg1
            if x2<0:x2=0
            ptg1=G1[2]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[2])-SDSL,int(G1[2])+SDSL)
                ptg1=random.choice(ptg)
            x3=ptg1
            if x3<0:x3=0
            ptg1=G1[3]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[3])-SDSL,int(G1[3])+SDSL)
                ptg1=random.choice(ptg)
            x4=ptg1
            if x4<0:x4=0
            ptg1=G1[4]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[4])-SDspeed,int(G1[4])+SDspeed)
                ptg1=random.choice(ptg)
            x5=ptg1
            if x5<0:x5=0
            ptg1=G1[5]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[5])-SDspeed,int(G1[5])+SDspeed)
                ptg1=random.choice(ptg)
            x6=ptg1
            if x6<0:x6=0
            ptg1=G1[6]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[6])-SDSL,int(G1[6])+SDSL)
                ptg1=random.choice(ptg)
            x7=ptg1
            if x7<0:x7=0
            ptg1=G1[7]
            mutation_rate2=random.choice(mutation_rate)
            if mutation_rate2=="y":
                ptg=range(int(G1[7])-SDSL,int(G1[7])+SDSL)
                ptg1=random.choice(ptg)
            x8=ptg1
            if x8<0:x8=0
            #taille de pistil
            size=x3+x4
            #nombre d'ovule
            OVUL2=round(OVUL-cost*size,0)
            if OVUL2<=0:OVUL2=0
            # nomdre d'ovule survivants
            survi=1
            # nombre de pollen acceuilli
            partner=ratio*OVUL
            G1=[x1,x2,x3,x4,x5,x6,x7,x8,size,OVUL2,survi,partner]
            creation1.append(G1)
            i=i+1
            
        if len(creation1)==0:
            generation_ini=nb_generation+1
            x="\n"+str(l)+";extinction"
            i=0
            while i <= len(creation1)-1:
                ind=creation1[i]
                fichier2.write(str(x)+";"+str(ind)) 
                i=i+1
            x=str(x)+";NA;NA;NA;NA;"+str(len(graine))+";NA;NA;NA;NA"
            fichier.write(str(x))
        
        if len(creation1)!=0 and generation_ini%10==0:        
            x="\n"+str(l)+";"+str(generation_ini)
            i=0
            SPEED=0
            STYLE=0
            spe_neu=0
            sty_neu=0
            while i <= len(creation1)-1:
                ind=creation1[i]
                fichier2.write(str(x)+";"+str(ind)) 
                SPEED=SPEED+ind[0]+ind[1]
                STYLE=STYLE+ind[2]+ind[3]
                spe_neu=spe_neu+ind[4]+ind[5]
                sty_neu=sty_neu+ind[6]+ind[7]
                i=i+1
            SPEED=SPEED/(2*N)
            STYLE=STYLE/(2*N)
            spe_neu=spe_neu/(2*N)
            sty_neu=sty_neu/(2*N)
            
            x=str(x)+";"+str(SPEED)+";"+str(STYLE)+";"+str(spe_neu)+";"+str(sty_neu)+";"+str(len(graine))
            fichier.write(str(x)) 
        suivi1=creation1[:]
        generation_ini=generation_ini+1
        
    print(l)
    l=l+1
            
fichier.close()
fichier2.close()
print("ok")

