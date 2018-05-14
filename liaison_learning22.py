#!/usr/bin/python
#Have eta proportionally less for h-aspire words
#In this version, we don't try to deal with optional enchainement.
#In this verson we model differing frequencies and their effect on liaison.

from __future__ import division
import random
import numpy as np
import json

#To start, ami has free vairants ami, tami, zami, nami
#joli starts out as joli

#petit starts out as peti

#Suppose that we start with faithfulness ranked low.
#Suppose that a learner encounters petit followed by a vowel initial word more often than ami encounters a preceding word with t as opposed to z or n. So we increment three times as often for petit. So we change t2 at a third the rate.


class WeightsAndActivations:

    def __init__(self):
        faithfactor = 0.5
        self.dep = 1.0*faithfactor
        self.maxx = 1.0*faithfactor
        self.unif = 1.0*faithfactor
        self.alignlC = 1.0*faithfactor
        self.alignrC = 1.0*faithfactor
        self.alignlV = 1.0*faithfactor
        self.alignrV = 1.0*faithfactor
        self.t2 = 0.3
        self.z2 =  0.3
        self.n2 =  0.3
        self.beta = 0.5
        self.beta2 = 0.5
        self.omega = 0.5
        self.omega2 = 0.5
        self.omega3 = 0.5
        self.phi = 0.5
        self.combined = [0.0, 0.0, 0.0]
        self.t1 = 0.3
        self.z1 = 0.3
        self.zs = 0.3
        self.n1 = 0.3
        self.t1f = 0.3
        self.z1f = 0.3
        self.n1f = 0.3
        self.null = 0
        self.ele = 1.0
        self.eha = 0.0
        self.eho = 0.0
        self.eham = 0.0
        self.ehe = 0.0
        self.onset = 1.0
        self.nocoda = 1.0
        self.complex = 1.0
        self.integrity = 1.0
        self.contiguity = 1.0
        self.mu = 1.0
        self.initialdep = 1.0*faithfactor
        self.initialmaxx = 1.0*faithfactor
        self.initialunif = 1.0*faithfactor
        self.initialalignlC = 1.0*faithfactor
        self.initialalignrC = 1.0*faithfactor
        self.initialalignlV = 1.0*faithfactor
        self.initialalignrV = 1.0*faithfactor
        self.initialt2 = 0.3
        self.initialz2 = 0.3
        self.initialn2 = 0.3
        self.initialbeta = 0.5
        self.initialbeta2 = 0.5
        self.initialomega = 0.5
        self.initialomega2 = 0.5
        self.initialomega3 = 0.5
        self.initialphi = 0.5
        self.initialcombined = [0.0, 0.0, 0.0]
        self.initialt1 = 0.3
        self.initialz1 = 0.3
        self.initialzs = 0.3
        self.initialn1 = 0.3
        self.initialt1f = 0.3
        self.initialz1f = 0.3
        self.initialn1f = 0.3
        self.initialnull = 0
        self.initialele = 1.0
        self.initialeha = 0.0
        self.initialeho = 0.0
        self.initialeham = 0.0
        self.initialehe = 0.0
        self.initialonset = 1.0
        self.initialnocoda = 1.0
        self.initialcomplex = 1.0
        self.initialintegrity = 1.0
        self.initialcontiguity = 1.0
        self.initialmu = 1.0
        return

    #Randomly choose one of three possible liaison consonants.
    def chooseL1(self):
        idx = random.randint(0,3)
        if idx == 0:
            return (self.null, 'null')
        elif idx == 1:
            return (self.t2, 't2')
        elif idx == 2:
            return (self.z2, 'z2')
        else:
            return (self.n2, 'n2')

    #Increment functions for various activations and weights:
    def inct2(self, increment):
        self.t2 += increment/3
        return

    def incz2(self, increment):
        self.z2 += increment/3
        return

    def incn2(self, increment):
        self.n2 += increment/3
        return

    def inccomb(self, increment, idx):
        self.combined[idx] += increment
        return

    def redcomb(self,factor):
        for i in range(3):
            self.combined[i] = self.combined[i]/factor

    def inct1(self, increment):
        self.t1 += increment
        return

    def incz1(self, increment):
        self.z1 += increment
        return

    def inczs(self, increment):
        self.zs += increment*7
        return

    def incn1(self, increment):
        self.n1 += increment
        return

    def inct1f(self, increment):
        self.t1f += increment
        return

    def incn1f(self, increment):
        self.n1f += increment
        return

    def incz1f(self, increment):
        self.z1f += increment
        return

    def incbeta(self, increment):
        self.beta += increment
        return

    def incbeta2(self, increment):
        self.beta2 += increment
        return

    def incomega(self, increment):
        self.omega += increment
        return

    def incomega2(self, increment):
        self.omega2 += increment
        return

    def incomega3(self, increment):
        self.omega3 += increment
        return

    def incphi(self, increment):
        self.phi += increment
        return

    def incele(self, increment):
        self.ele += increment
        return

    def inceha(self, increment):
        self.eha += increment
        return

    def inceho(self, increment):
        self.eho += increment
        return

    def inceham(self, increment):
        self.eham += increment
        return

    def incehe(self, increment):
        self.ehe += increment
        return

    def incdep(self, increment):
        self.dep += increment
        return

    def incmax(self, increment):
        self.maxx += increment
        return

    def inconset(self, increment):
        self.onset += increment
        return

    def incunif(self, increment):
        self.unif += increment
        return

    def inccomplex(self, increment):
        self.complex += increment
        return

    def inccontig(self, increment):
        self.contiguity += increment
        return

    def incalignlC(self, increment):
        self.alignlC += increment/2
        return

    def incalignrC(self, increment):
        self.alignrC += increment/10
        return

    def incalignlV(self, increment):
        self.alignlV += increment
        return

    def incalignrV(self, increment):
        self.alignrV += increment
        return

    def incnocoda(self, increment):
        self.nocoda += increment/2
        return

    def decrmu(self, iteration_number):
        self.mu = 1.0 / (1.0 + iteration_number)
        return

    def setehazero(self):
        self.eha = 0.0
        return

    def setehezero(self):
        self.ehe = 0.0
        return

wa = WeightsAndActivations()

#Stepsizes for learing:
eta0 = 0.03
eta1 = 0.03
eta2 = 0.05
margin = 0.03
boundary_limit = 0.9
combined_limit = 0.05

'''What happens if the learner chooses the wrong allomorph -- e.g. un tami  or petit zami. In that case, no candidate will match the correct form that the child hears. Suppose then that the child adds the correct consonant to the adjective.'''

for i in range(500):
    all_correct = True
    print 'iteration', i
    print

    #petit vs. peti They both satisfy AlignrV so we leave it out.
    h_peti = 0
    h_petit = wa.alignrC + wa.maxx*min(1.0,wa.t1) - wa.dep*(1.0 - min(wa.t1,1)) -wa.nocoda
    print 'Harmonies for petit.'
    print 'peti', h_peti, 'petit', h_petit
    if h_peti - margin > h_petit:
        print 'Correct candidate for peti'
    else:
        all_correct = False
        print 'Wrong winner for peti. Winner should be the first value' 
        #Decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease t1
        if wa.t1 > eta2:
            wa.inct1(-eta2*3/4) #Because we have 4 examples that decrease t1 as opposed to 3 for the others
            print 'Decreased t1 by', eta2
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx > eta1:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
            

    print
    print


    #petite
    h_peti = 0
    h_petite = wa.alignrC + wa.maxx*min(1.0,wa.t1f) - wa.dep*(1.0 - min(wa.t1f,1)) - wa.nocoda

    print 'harmonies for petite'
    print 'petit', h_petite, 'peti', h_peti
    if h_petite -margin > h_peti:
        print 'Correct candidate for petite'
    else:
        all_correct = False
        print 'Wrong winner for petite. Winner should be the first value' 
        #increase alignr
        wa.incalignrC(eta1)
        #increase t1f
        if wa.t1f < 1-eta2:
            wa.inct1f(eta2)
        #decrease dep
        if wa.dep > eta1:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)

    print
    print
   
    #bonne
    h_bo = 0
    h_bonne = wa.alignrC + wa.maxx*min(1.0,wa.n1f) - wa.dep*(1.0 - min(wa.t1f,1)) - wa.nocoda
    print 'harmonies for bonne'
    print 'bon', h_bonne, 'bo', h_bo
    if h_bonne - margin > h_bo:
        print 'Correct candidate for bonne'
    else:
        all_correct = False
        print 'Wrong winner for bonne. Winner should be the first value' 
        #increase alignr
        wa.incalignrC(eta1)
        #increase n1f
        if wa.n1f < 1-eta2:
            wa.incn1f(eta2)
        #decrease dep
        if wa.dep > eta1:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)
           

    print
    print

   
    #douze
    h_dou = 0
    h_douze  = wa.alignrC + wa.maxx*min(1.0,wa.z1f) - wa.dep*(1.0 - min(wa.z1f,1)) - wa.nocoda
    print 'harmonies for douze'
    print 'douze', h_douze, 'dou', h_dou
    if h_douze - margin > h_dou:
        print 'Correct candidate for douze'
    else:
        all_correct = False
        print 'Wrong winner for douze. Winner should be the first value' 
        #increase alignr
        wa.incalignrC(eta1)
        #increase z1f
        if wa.z1f < 1-eta2:
            wa.incz1f(eta2)
        #decrease dep
        if wa.dep > eta1:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)
           

    print
    print
    #deux vs. deu
    h_deu = 0
    h_deuz = wa.alignrC + wa.maxx*min(1.0,wa.z1) - wa.dep*(1.0 - min(wa.z1,1)) -wa.nocoda
    print 'harmonies for deux'
    print 'deu', h_deu, 'deuz', h_deuz
    if h_deu -margin > h_deuz:
        print 'Correct candidate for deu'
    else:
        all_correct = False
        print 'Wrong winner for deu. Winner should be the first value' 
        #Decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease z1
        if wa.z1 > eta2:
            wa.incz1(-eta2)
            print 'Decreased z1 by', eta2
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx > eta1:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
            

    print
    print



    #bo vs. bon
    h_bo = 0
    h_bon = wa.alignrC + wa.maxx*min(1.0,wa.n1) - wa.dep*(1.0 - min(wa.n1,1)) -wa.nocoda
    print 'harmonies for bon'
    print 'bo', h_bo, 'bon', h_bon
    if h_bo - margin > h_bon:
        print 'Correct candidate for bo'
    else:
        all_correct = False
        print 'Wrong winner for bo. Winner should be the first value' 
        #Decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease n1
        if wa.n1 > eta2:
            wa.incn1(-eta2)
            print 'Decreased n1 by', eta2
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx > eta1:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
            
    print
    print

    #Suppose that the probability that they consider there to be three allomorphs as opposed to a blend is 0.8 to the power of the number of the iteration.
    prob = 0.8**i
    r = random.random()
    if r < prob:
        mult_allo = True
        (allomorph, allomorph_name) = wa.chooseL1()
    else:
        mult_allo = False
    
    #ami: input is one of tami, zami, nami, ami or else a blend They all satisfy AlignlV assuming no schwa activation on ami so we leave out AlignlV
    print 'Checking ami'
    h_ami = -wa.onset
    if mult_allo:
        if allomorph_name == 't2':
            print 'Allomorph is t2.'
            h_t2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_z2ami = -100.0
            h_n2ami = -100.0
        elif allomorph_name == 'z2':
            print 'Allomorph is z2.'
            h_z2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_t2ami = -100.0
            h_n2ami = -100.0
        elif allomorph_name == 'n2':
            print 'Allomorph is n2.'
            h_n2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_t2ami = -100.0
            h_z2ami = -100.0
        else:
            print 'Allomorph is null.'
            h_t2ami = -100.0
            h_z2ami = -100.0
            h_n2ami = -100.0
    else:
        print 'Blended t2, z2, n2'
        h_t2ami =  -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC
        h_z2ami =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC
        h_n2ami =  -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC
    print 'harmonies for ami'
    print 'ami', h_ami, 'tami', h_t2ami, 'zami', h_z2ami, 'nami', h_n2ami
    harmonies = [h_ami - margin, h_t2ami, h_z2ami, h_n2ami]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be ami
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for ami'
    else:
        print 'Incorrect winner for ami. Highest is ', highest
        all_correct = False
        if 1 in highest:
            #reduce activation on t2
            if wa.t2 > eta2:
                wa.inct2(-eta2)
        if 2 in highest:
            #reduce activation on z2
            if wa.z2 > eta2:
                wa.incz2(-eta2)
        if 3 in highest:
            #reduce activation on n2
            if wa.n2 > eta2:
                wa.incn2(-eta2)
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #decrease onset
        if wa.onset > eta1+0.01:
            wa.inconset(-eta1)
        #decrease alignl
        if wa.alignlC > eta1:
            wa.incalignlC(-eta1)
            print 'Decreasing alignl.'


    print
    print

    print 'Checking petit ami'
    h_petiami = -wa.onset
    #t1 may start out as zero.
    h_petit1ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1)
    h_petit1_ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that z1 and n1 will not occur with petit
    h_petizami = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_petiz_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    h_petinami = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_petin_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 't2':
            #Input is peti(t1) + t2ami
            print 'Allomorph is t2.'
            h_petit2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_petit2_ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_petit12ami = -wa.dep*(1.0 - min((wa.t1 + allomorph),1)) + wa.maxx*min(1,(wa.t1 + allomorph)) + wa.alignlC - wa.unif*wa.mu
            h_petit12_ami = -wa.dep*(1.0 - min((wa.t1 + allomorph),1)) + wa.maxx*min(1,(wa.t1 + allomorph)) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
            h_petit1_z2ami = -100.0
            h_petit1_n2ami = -100.0
            h_petiz2ami = -100.0
            h_petiz2_ami = -100.0
            h_petin2ami = -100.0
            h_petin2_ami = -100.0
            #Are the next 4 even candidates if z1 and z2  n1 and n2 not there in input?
            #h_petiz12ami = -wa.dep - wa.unif*wa.mu #Not left anchored because z2 not there in input
            #h_petiz12_ami = -wa.dep + wa.alignl - wa.unif*wa.mu - wa.nocoda - wa.onset 
            #h_petin12ami = -wa.dep - wa.unif*wa.mu
            #h_petin12_ami = -wa.dep* + wa.alignl - wa.unif*wa.mu - wa.nocoda - wa.onset
        elif allomorph_name == 'z2':
            #input is peti(t1) + z2ami
            print 'Allomorph is z2.'
            h_petit2ami = -100.0
            h_petit2_ami = -100.0
            h_petit12ami = -100.0
            h_petit12_ami = -100.0
            h_petit1_z2ami =  -wa.dep*(1.0-min(wa.t1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t1) + wa.maxx*min(1,wa.z2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_petit1_n2ami = -100.0
            h_petiz2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  + wa.alignlC 
            h_petiz2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  - wa.nocoda - wa.onset
            h_petin2ami = -100.0
            h_petin2_ami = -100.0
            #There should be no colaseced candidates
            #h_petiz12ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignl - wa.unif*wa.mu
            #h_petiz12_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignr - wa.unif*wa.mu - wa.nocoda - wa.onset
            #h_petit12ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1) - wa.unif*wa.mu
            #h_petit12_ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1) -wa.onset - wa.nocoda -wa.unif*wa.mu
            #h_petin12ami =  -wa.dep + wa.alignr - wa.unif*wa.mu #This is anchored right if n1 and n2 are not in the input
        elif allomorph_name == 'n2':
            #input is peti(t1) + n2ami
            print 'Allomorph is n2.'
            h_petit2ami = -100.0
            h_petit2_ami = -100.0
            h_petit12ami = -100.0
            h_petit12_ami = -100.0
            h_petit1_z2ami = -100.0
            h_petit1_n2ami =  -wa.dep*(1.0-min(wa.t1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t1) + wa.maxx*min(1,wa.n2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_petiz2ami = -100.0
            h_petiz2_ami = -100.0
            h_petin2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC 
            h_petin2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
            #There should be no colaseced candidates
            #h_petit12ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1)
            #h_petin12ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignl - wa.unif*wa.mu
            #h_petiz12ami =  -wa.dep + wa.alignl - wa.unif*wa.mu
        else:
            #input is petit(t1) + ami
            print 'Allomorph is null.'
            h_petit2ami = -100.0
            h_petit2_ami = -100.0
            h_petit12ami = -100.0
            h_petit12_ami = -100.0
            h_petit1_z2ami = -100.0
            h_petit1_n2ami = -100.0
            h_petiz2ami = -100.0
            h_petiz2_ami = -100.0
            h_petin2ami = -100.0
            h_petin2_ami = -100.0
            #There should be no colaseced candidates
            #h_petit12ami = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1)
            #h_petiz12ami = -wa.dep + wa.alignl - wa.unif*wa.mu
            #h_petin12ami = -wa.dep + wa.alignl - wa.unif*wa.mu
    else:
        print 'Blended t2, z2, n2'
        #input is peti(t1) + [t2,z2,n2]ami
        h_petit2ami =  -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC
        h_petit2_ami =  -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_petit12ami = -wa.dep*(1.0 - min((wa.t1 + wa.combined[0]),1)) + wa.maxx*min(1,(wa.t1 + wa.combined[0])) + wa.alignlC - wa.unif*wa.mu
        h_petit12_ami = -wa.dep*(1.0 - min((wa.t1 + wa.combined[0]),1)) + wa.maxx*min(1,(wa.t1 + wa.combined[0])) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
        h_petit1_z2ami = -100.0
        h_petit1_n2ami = -100.0
        #These next two should not be candidates if there is no z1 or n1 in the input
        #h_petiz12ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignl - wa.unif*wa.mu
        #h_petin12ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignl - wa.unif*wa.mu
        h_petiz2ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_petiz2_ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_petin2ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_petin2_ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset

    print 'harmonies for petitami'
    print h_petiami, h_petit1ami, h_petit1_ami, h_petizami, h_petiz_ami, h_petinami, h_petin_ami, h_petit2ami, h_petit2_ami, h_petit12ami, h_petit12_ami, h_petit1_z2ami, h_petit1_n2ami, h_petiz2ami, h_petin2ami, h_petiz2_ami, h_petin2_ami, h_petit1_z2ami 

#                0            1            2           3            4

    harmonies = [ h_petiami, h_petit1ami, h_petit1_ami, h_petizami, h_petiz_ami, h_petinami, h_petin_ami, h_petit2ami, h_petit2_ami, h_petit12ami, h_petit12_ami, h_petit1_z2ami, h_petit1_n2ami, h_petiz2ami, h_petin2ami, h_petiz2_ami, h_petin2_ami, h_petit1_z2ami, h_petit1_z2ami] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible members of highest to be petit1ami and petit12ami and petit2ami and petit1_ami and petit12_ami and petit2_ami
    print 'dep:', wa.dep, 't1:', wa.t1, 't2', wa.t2, 'combined0', wa.combined[0], 'max', wa.maxx, 'alignlC', wa.alignlC, 'unif', wa.unif*wa.mu 
    bad = [0,3,4,5,6,11,12,13,14,15,16,17,18]
    if list(set(highest) & set(bad)) == []:
    #if any of 1,2,7,8,9,10 in highest:
        print 'For petitami, highest exclusively includes at least one correct candidate and should be any of 1,2,7,8,9,10', highest
    else:
        all_correct = False
        print 'Wrong result for petitami. Winning candidate indices should be exclusive from the set 1,2,u7,8,9,10 and are are', highest
        #We want petit1ami or petit12ami to win. We can do the following.
        #increase t1
        wa.inct1(eta2/2)
        print 'increased t1'
        #increase t2
        if mult_allo:
            wa.inct2(eta2/2)
            print 'increased t2'
        else:
            wa.inccomb(eta2/2,0)
            print 'increased comb0'
        #decrease uniformity
        if wa.unif*wa.mu > eta0:
            wa.incunif(-eta0)
        #Increase alignl
        wa.incalignlC(eta1)
        print 'increasing alignl.'

        #decrease dep as long as it doesn't go below zero
        if 0 in highest: 
            #decrease dep
            if wa.dep > eta1*5/13:
                wa.incdep(-eta1*5/13)
            #increase onset
            wa.inconset(eta1)
            #increase maxx
            wa.incmax(eta1*5/13)
            print 'decreasing dep and increasing max'
        if list(set(highest) & set([11,13,15,17,18])) != []:
            #decrease z2
            if mult_allo:
                if wa.z2 >= eta2/2:
                    wa.incz2(-eta2/2)
            else:
                if wa.combined[1] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,1)
                    print 'decreased comb1', wa.combined[1]
        if list(set(highest) & set([12,14,16,18])) != []:
            #decrease n2
            if mult_allo:
                if wa.n2 >= eta2/2:
                    wa.incn2(-eta2/2)
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]

        if list(set(highest) & set([2,4,6,8,10,11,12,15,16,17,18])) != []:
            wa.incnocoda(eta1)
            wa.inconset(eta1)

    print



    print
    print

    print 'Checking bon ami.'
    h_boami = -wa.onset
    #t1 may start out as zero.
    h_bon1ami = -wa.dep*(1.0 - min(wa.n1,1)) + wa.maxx*min(1,wa.n1)
    h_bon1_ami = -wa.dep*(1.0 - min(wa.n1,1)) + wa.maxx*min(1,wa.n1) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that z1 and t1 will not occur with bon
    h_botami = -wa.dep #Not right anchored as long as n1 is there and also not left anchored
    h_bot_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as n1 is there 
    h_bozami = -wa.dep #Not right anchored as long as n1 is there and also not left anchored
    h_boz_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as n1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'n2':
            #Input is bo(n1) + n2ami
            print 'Allomorph is n2.'
            h_bon2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_bon2_ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_bon12ami = -wa.dep*(1.0 - min((wa.n1 + allomorph),1)) + wa.maxx*min(1,(wa.n1 + allomorph)) + wa.alignlC - wa.unif*wa.mu
            h_bon12_ami = -wa.dep*(1.0 - min((wa.n1 + allomorph),1)) + wa.maxx*min(1,(wa.n1 + allomorph)) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
            h_bon1_t2ami = -100.0
            h_bon1_z2ami = -100.0
            h_bot2ami = -100.0
            h_bot2_ami = -100.0
            h_boz2ami = -100.0
            h_boz2_ami = -100.0
        elif allomorph_name == 'z2':
            #input is bo(n1) + z2ami
            print 'Allomorph is z2.'
            h_bon2ami = -100.0
            h_bon2_ami = -100.0
            h_bon12ami = -100.0
            h_bon12_ami = -100.0
            h_bon1_z2ami =  -wa.dep*(1.0-min(wa.n1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n1) + wa.maxx*min(1,wa.z2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_bon1_t2ami = -100.0
            h_boz2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  + wa.alignlC 
            h_boz2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  - wa.nocoda - wa.onset
            h_bot2ami = -100.0
            h_bot2_ami = -100.0
        elif allomorph_name == 't2':
            #input is bo(n1) + t2ami
            print 'Allomorph is t2.'
            h_bon2ami = -100.0
            h_bon2_ami = -100.0
            h_bon12ami = -100.0
            h_bon12_ami = -100.0
            h_bon1_z2ami = -100.0
            h_bon1_t2ami =  -wa.dep*(1.0-min(wa.n1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n1) + wa.maxx*min(1,wa.t2) + wa.alignrC - wa.alignlC - wa.nocoda 
            h_boz2ami = -100.0
            h_boz2_ami = -100.0
            h_bot2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC 
            h_bot2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset
        else:
            #input is bon(n1) + ami
            print 'Allomorph is null.'
            h_bon2ami = -100.0
            h_bon2_ami = -100.0
            h_bon12ami = -100.0
            h_bon12_ami = -100.0
            h_bon1_z2ami = -100.0
            h_bon1_t2ami = -100.0
            h_bot2ami = -100.0
            h_bot2_ami = -100.0
            h_boz2ami = -100.0
            h_boz2_ami = -100.0
    else:
        print 'Blended t2, z2, n2'
        #input is bo(n1) + [t2,z2,n2]ami
        h_bon2ami =  -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC
        h_bon2_ami =  -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset
        h_bon12ami = -wa.dep*(1.0 - min((wa.n1 + wa.combined[2]),1)) + wa.maxx*min(1,(wa.n1 + wa.combined[2])) + wa.alignlC - wa.unif*wa.mu
        h_bon12_ami = -wa.dep*(1.0 - min((wa.n1 + wa.combined[2]),1)) + wa.maxx*min(1,(wa.n1 + wa.combined[2])) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
        h_bon1_z2ami = -100.0
        h_bon1_t2ami = -100.0
        #These next two should not be candidates if there is no z1 or n1 in the input
        #h_bot12ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignl - wa.unif*wa.mu
        #h_boz12ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignl - wa.unif*wa.mu
        h_bot2ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_bot2_ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_boz2ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_boz2_ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset

    print 'harmonies for bon ami'
    print h_boami, h_bon1ami, h_bon1_ami, h_botami, h_bot_ami, h_bozami, h_boz_ami, h_bon2ami, h_bon2_ami, h_bon12ami, h_bon12_ami, h_bon1_z2ami, h_bon1_t2ami, h_bot2ami, h_boz2ami, h_bot2_ami, h_boz2_ami

#                0            1           2         3         4

    harmonies =[ h_boami, h_bon1ami, h_bon1_ami, h_botami, h_bot_ami, h_bozami, h_boz_ami, h_bon2ami, h_bon2_ami, h_bon12ami, h_bon12_ami, h_bon1_z2ami, h_bon1_t2ami, h_boz2ami, h_bot2ami, h_boz2_ami, h_bot2_ami]  #The last is a duplicate
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible members of highest to be bon1ami and bon12ami and bon2ami and bon1_ami and bon12_ami and bon2_ami
    bad = [0,3,4,5,6,11,12,13,14,15,16,17,18]
    if list(set(highest) & set(bad)) == []:
    #if any of 1,2,7,8,9,10 in highest:
        print 'For bonami, highest exclusively includes at least one correct candidate and should be any of 1,7,9', highest
    else:
        all_correct = False
        print 'Wrong result for bonami. Winning candidate indices should be exclusive from the set 1,7,9 and are are', highest
        #We want bon1ami or bon12ami to win. We can do the following.
        #increase n1
        wa.incn1(eta2/2)
        print 'increased n1'
        #increase n2
        if mult_allo:
            wa.incn2(eta2/2)
            print 'increased n2'
        else:
            wa.inccomb(eta2/2,2)
            print 'increased comb2'
        #decrease uniformity
        if wa.unif*wa.mu > eta0:
            wa.incunif(-eta0)
        #Increase alignl
        wa.incalignlC(eta1)
        print 'increasing alignl.'

        #decrease dep as long as it doesn't go below zero
        if 0 in highest: 
            #decrease dep
            if wa.dep > eta1*5/13:
                wa.incdep(-eta1*5/13)
            #increase onset
            wa.inconset(eta1)
            #increase maxx
            wa.incmax(eta1*5/13)
            print 'decreasing dep and increasing max'
            #decrease uniformity
            #if wa.unif*wa.mu > eta0:
             #   wa.incunif(-eta0)
        if list(set(highest) & set([11,13,15,17])) != []:
            #decrease z2
            if mult_allo:
                if wa.z2 >= eta2/2:
                    wa.incz2(-eta2/2)
                    print 'decreased z2'
            else:
                if wa.combined[1] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,1)
                    print 'decreased comb1', wa.combined[1]
        if list(set(highest) & set([12,14,16,18])) != []:
            #decrease t2
            if mult_allo:
                if wa.t2 >= eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2'
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased comb0', wa.combined[0]

        if list(set(highest) & set([2,4,6,8,10,11,12,15,16])) != []:
            wa.incnocoda(eta1)
            wa.inconset(eta1)

    print




    print
    print

    print 'Checking deux amis'
    h_deuami = -wa.onset
    #t1 may start out as zero.
    h_deuz1ami = -wa.dep*(1.0 - min(wa.z1,1)) + wa.maxx*min(1,wa.z1)
    h_deuz1_ami = -wa.dep*(1.0 - min(wa.z1,1)) + wa.maxx*min(1,wa.z1) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that n1 and t1 will not occur with deuz
    h_deunami = -wa.dep #Not right anchored as long as z1 is there and also not left anchored
    h_deun_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as z1 is there 
    h_deutami = -wa.dep #Not right anchored as long as z1 is there and also not left anchored
    h_deut_ami = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as z1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'z2':
            #Input is deu(z1) + z2ami
            print 'Allomorph is z2.'
            h_deuz2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_deuz2_ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_deuz12ami = -wa.dep*(1.0 - min((wa.z1 + allomorph),1)) + wa.maxx*min(1,(wa.z1 + allomorph)) + wa.alignlC - wa.unif*wa.mu
            h_deuz12_ami = -wa.dep*(1.0 - min((wa.z1 + allomorph),1)) + wa.maxx*min(1,(wa.z1 + allomorph)) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
            h_deuz1_t2ami = -100.0
            h_deuz1_n2ami = -100.0
            h_deun2ami = -100.0
            h_deun2_ami = -100.0
            h_deut2ami = -100.0
            h_deut2_ami = -100.0
        elif allomorph_name == 't2':
            #input is deu(z1) + t2ami
            print 'Allomorph is t2.'
            h_deuz2ami = -100.0
            h_deuz2_ami = -100.0
            h_deuz12ami = -100.0
            h_deuz12_ami = -100.0
            h_deuz1_t2ami =  -wa.dep*(1.0-min(wa.z1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z1) + wa.maxx*min(1,wa.t2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_deuz1_n2ami = -100.0
            h_deut2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC 
            h_deut2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset
            h_deun2ami = -100.0
            h_deun2_ami = -100.0
        elif allomorph_name == 'n2':
            #input is deu(z1) + n2ami
            print 'Allomorph is n2.'
            h_deuz2ami = -100.0
            h_deuz2_ami = -100.0
            h_deuz12ami = -100.0
            h_deuz12_ami = -100.0
            h_deuz1_t2ami = -100.0
            h_deuz1_n2ami =  -wa.dep*(1.0-min(wa.z1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z1) + wa.maxx*min(1,wa.n2) + wa.alignrC - wa.alignlC - wa.nocoda 
            h_deut2ami = -100.0
            h_deut2_ami = -100.0
            h_deun2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC 
            h_deun2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
        else:
            #input is deuz(z1) + ami
            print 'Allomorph is null.'
            h_deuz2ami = -100.0
            h_deuz2_ami = -100.0
            h_deuz12ami = -100.0
            h_deuz12_ami = -100.0
            h_deuz1_n2ami = -100.0
            h_deuz1_t2ami = -100.0
            h_deun2ami = -100.0
            h_deun2_ami = -100.0
            h_deut2ami = -100.0
            h_deut2_ami = -100.0
    else:
        print 'Blended t2, z2, n2'
        #input is deu(z1) + [t2,z2,n2]ami
        h_deuz2ami =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC
        h_deuz2_ami =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_deuz12ami = -wa.dep*(1.0 - min((wa.z1 + wa.combined[1]),1)) + wa.maxx*min(1,(wa.z1 + wa.combined[1])) + wa.alignlC - wa.unif*wa.mu
        h_deuz12_ami = -wa.dep*(1.0 - min((wa.z1 + wa.combined[1]),1)) + wa.maxx*min(1,(wa.z1 + wa.combined[1])) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
        h_deuz1_n2ami = -100.0
        h_deuz1_t2ami = -100.0
        h_deun2ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_deun2_ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset
        h_deut2ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_deut2_ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset

    print 'harmonies for deux amis'
    print h_deuami, h_deuz1ami, h_deuz1_ami, h_deunami, h_deun_ami, h_deutami, h_deut_ami, h_deuz2ami, h_deuz2_ami, h_deuz12ami, h_deuz12_ami, h_deuz1_n2ami, h_deuz1_t2ami, h_deun2ami, h_deut2ami, h_deun2_ami, h_deut2_ami

#                0            1           2            3           4

    harmonies =[ h_deuami, h_deuz1ami, h_deuz1_ami, h_deunami, h_deun_ami, h_deutami, h_deut_ami, h_deuz2ami, h_deuz2_ami, h_deuz12ami, h_deuz12_ami, h_deuz1_n2ami, h_deuz1_t2ami, h_deun2ami, h_deut2ami, h_deun2_ami, h_deut2_ami] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible members of highest to be deuz1ami and deuz12ami and deuz2ami and deuz1_ami and deuz12_ami and deuz2_ami
    bad = [0,3,4,5,6,11,12,13,14,15,16,17,18]
    if list(set(highest) & set(bad)) == []:
    #if any of 1,2,7,8,9,10 in highest:
        print 'For deuzami, highest exclusively includes at least one correct candidate and should be any of 1,7,9', highest
    else:
        all_correct = False
        print 'Wrong result for deuzami. Winning candidate indices should be exclusive from the set 1,7,90 and are are', highest
        #We want deuz1ami or deuz12ami to win. We can do the following.
        #increase z1
        wa.incz1(eta2/2)
        #increase z2
        if mult_allo:
            wa.incz2(eta2/2)
            print 'increased z2'
        else:
            wa.inccomb(eta2/2,1)
            print 'increased comb1'
        #decrease uniformity
        if wa.unif*wa.mu > eta0:
            wa.incunif(-eta0)
        #Increase alignl
        wa.incalignlC(eta1)
        print 'increasing alignl.'

        #decrease dep as long as it doesn't go below zero
        if 0 in highest: 
            #decrease dep
            if wa.dep > eta1*5/13:
                wa.incdep(-eta1*5/13)
            #increase onset
            wa.inconset(eta1*5/13)
            #increase maxx
            wa.incmax(eta1)
            print 'decreasing dep and increasing max'
        if list(set(highest) & set([11,13,15])) != []:
            #decrease n2
            if mult_allo:
                if wa.n2 >= eta2/2:
                    wa.incn2(-eta2/2)
                    print 'decreased n2'
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]
        if list(set(highest) & set([12,14,16])) != []:
            #decrease t2
            if mult_allo:
                if wa.t2 >= eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2'
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased comb0', wa.combined[0]

        if list(set(highest) & set([2,4,6,8,10,11,12,15,16,17,18])) != []:
            wa.incnocoda(eta1)
            wa.inconset(eta1)
        #if list(set(highest) & set([2,10,11,12])) != []:
        #    if wa.alignrC > eta1:
        #        wa.incalignrC(-eta1)

    print

    print 'Checking joli ami'
    h_joliami = -wa.onset + wa.alignrV
    #Assume that we never give activation to a final consonant on joli
    #Assume that t1 and z1 and n1 will not occur with joli
    #THe candidates with VC.V syllabilification can't possibly win but are put in for completeness. No learning is incldued for these as false winners since they can't be.
    if mult_allo:
        (allomorph, allomorph_name) = wa.chooseL1()
        if allomorph_name == 't2':
            print 'Allomorph is t2.'
            h_jolit2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrV
            h_joliz2ami = -wa.dep + wa.alignrV
            h_jolin2ami = -wa.dep + wa.alignrV
            h_jolit2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda
            h_joliz2_ami = -wa.dep -wa.onset - wa.nocoda
            h_jolin2_ami = -wa.dep -wa.onset - wa.nocoda
        elif allomorph_name == 'z2':
            print 'Allomorph is z2.'
            h_joliz2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrV
            h_jolit2ami = -wa.dep + wa.alignrV
            h_jolin2ami = -wa.dep + wa.alignrV
            h_joliz2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda
            h_jolit2_ami = -wa.dep -wa.onset - wa.nocoda
            h_jolin2_ami = -wa.dep -wa.onset - wa.nocoda
        elif allomorph_name == 'n2':
            print 'Allomorph is n2.'
            h_jolin2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrV
            h_jolit2ami = -wa.dep + wa.alignrV
            h_joliz2ami = -wa.dep + wa.alignrV
            h_jolin2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda
            h_jolit2_ami = -wa.dep -wa.onset - wa.nocoda
            h_joliz2_ami = -wa.dep -wa.onset - wa.nocoda
        else:
            print 'Allomorph is null.'
            h_jolit2ami = -wa.dep + wa.alignrV
            h_joliz2ami = -wa.dep + wa.alignrV
            h_jolin2ami = -wa.dep + wa.alignrV
            h_jolit2_ami = -wa.dep -wa.onset - wa.nocoda
            h_joliz2_ami = -wa.dep -wa.onset - wa.nocoda
            h_jolin2_ami = -wa.dep -wa.onset - wa.nocoda
    else:
        h_jolit2ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC + wa.alignrV
        h_joliz2ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC + wa.alignrV
        h_jolin2ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC + wa.alignrV
        h_jolit2_ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) -wa.onset - wa.nocoda
        h_joliz2_ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) -wa.onset - wa.nocoda
        h_jolin2_ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) -wa.onset - wa.nocoda


    print 'Harmonies for joli ami.'
    print h_joliami, h_jolit2ami, h_joliz2ami, h_jolin2ami, h_jolit2_ami, h_joliz2_ami, h_jolin2_ami
    harmonies = [h_joliami-margin, h_jolit2ami, h_joliz2ami, h_jolin2ami, h_jolit2_ami, h_joliz2_ami, h_jolin2_ami]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be joliami
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for joliami'
    else:
        all_correct = False
        print 'Wrong result for joliami. Winning index should be 0 and are:', highest
        if 0 not in highest:
            #Decrease onset.
            if wa.onset >= eta1:
                wa.inconset(-eta1)
            #increase dep
            wa.incdep(eta1*5/13)
            #decrease max
            if wa.maxx > eta1:
                wa.incmax(-eta1*5/13)
            print 'increasing dep and decreasing max'
            #decrease alignl
            if wa.alignlC > eta1:
                wa.incalignlC(-eta1)
                print 'decreasing alignl.'

        if 1 in highest:
            #decrease t2
            if mult_allo:
                if wa.t2 > eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2 by', eta2/2
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased combined[0] by', eta2/2, wa.combined[0]

        if 2 in highest:
            #decrease z2
            if mult_allo:
                if wa.z2 > eta2/2:
                    wa.incz2(-eta2/2)
                    print 'decreased z2'
            else:
                if wa.combined[1] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,1)
                    print 'decreased comb1', wa.combined[1]

        if 3 in highest:
            #decrease n2
            if mult_allo:
                if wa.n2 > eta2/2:
                    wa.incn2(-eta2/2)
                    print 'decreased n2'
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]

        if list(set(highest) & set([4,5,6])) != []:
            wa.incnocoda(eta1)


    print

    print
    print
    print
    #petite copine
    #input is peti(t1f)+copine

    print 'Checking petite copine.'
    h_peticopine = wa.alignlC
    h_petit1fcopine = -wa.dep*(1-wa.t1f) + wa.maxx*wa.t1f - wa.nocoda + wa.alignlC + wa.alignrC
    print h_petit1fcopine, h_peticopine
    harmonies = [h_petit1fcopine - margin, h_peticopine] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be petitecopine
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for petitecopine'
    else:
        all_correct = False
        print 'Incorrect result for petitecopine. Winning candidate should be 0 and indices are', highest
        #decrease dep
        if wa.dep >eta1*5/13:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)
        #increase alignr
        wa.incalignrC(eta1)
        #increase t1f
        if wa.t1f < 1-eta2:
            wa.inct1f(eta2)

    print
    print


    print
    print 'Checking deux copains.'
    #deux copains
    h_deucopain = wa.alignlC
    h_deuz1copain = -wa.dep*(1-wa.z1) + wa.maxx*wa.z1 - wa.nocoda + wa.alignlC + wa.alignrC
    print h_deucopain, h_deuz1copain 
    harmonies = [h_deucopain-margin, h_deuz1copain] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be deucopain
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for deucopain'
    else:
        all_correct = False
        print 'Incorrect result for deucopain. Winning candidate should be 0 and indices are', highest
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx >= eta1:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
        #decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease z1
        if wa.z1 > eta2:
            wa.incz1(-eta2)
            print 'decreased z1'

    print
    print


    print
    print 'Checking douze copains.'
    #douze copains
    h_doucopain = wa.alignlC
    h_douz1fcopain = -wa.dep*(1-wa.z1f) + wa.maxx*wa.z1f - wa.nocoda + wa.alignlC + wa.alignrC
    print h_douz1fcopain, h_doucopain 
    harmonies = [h_douz1fcopain - margin, h_doucopain] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be doucopain
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for douzecopain'
    else:
        all_correct = False
        print 'Incorrect result for douzecopain. Winning candidate should be 0 and indices are', highest
        #decrease dep
        if wa.dep >eta1*5/13:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)
        #increase alignr
        wa.incalignrC(eta1)
        #increase z1f
        if wa.z1f < 1-eta2:
            wa.incz1f(eta2)

    print
    print


    print
    print 'Checking bon copain.'
    #bon copain
    h_bocopain = wa.alignlC
    h_bon1copain = -wa.dep*(1-wa.n1) + wa.maxx*wa.n1 - wa.nocoda + wa.alignlC + wa.alignrC
    print h_bocopain, h_bon1copain 
    harmonies = [h_bocopain-margin, h_bon1copain] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be bocopain
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for bocopain'
    else:
        all_correct = False
        print 'Incorrect result for bocopain. Winning candidate should be 0 and indices are', highest
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx >= eta1:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
        #decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease n1
        if wa.n1 > eta2:
            wa.incn1(-eta2)
            print 'decreased n1'

    print
    print


    print
    print 'Checking bonne copine.'
    #bonne copine
    h_bocopine = wa.alignlC
    h_bon1fcopine = -wa.dep*(1-wa.n1f) + wa.maxx*wa.n1f - wa.nocoda + wa.alignlC + wa.alignrC
    print h_bon1fcopine, h_bocopine 
    harmonies = [h_bon1fcopine - margin, h_bocopine] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be bocopine
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for bonnecopine'
    else:
        all_correct = False
        print 'Incorrect result for bonnecopine. Winning candidate should be 0 and indices are', highest
        #decrease dep
        if wa.dep >eta1*5/13:
            wa.incdep(-eta1*5/13)
        #increase max
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #decrease nocoda
        if wa.nocoda > eta1:
            wa.incnocoda(-eta1)
        #increase alignr
        wa.incalignrC(eta1)
        #increase n1f
        if wa.n1f < 1-eta2:
            wa.incn1f(eta2)

    print
    print



    print
    print

    print 'Checking petit copain.'
    h_peticopain = wa.alignlC
    h_petit1copain = -wa.dep*(1-wa.t1) + wa.maxx*wa.t1 - wa.nocoda + wa.alignlC + wa.alignrC
    print h_peticopain, h_petit1copain  #, h_petiz1copain, h_petin1copain 
    harmonies = [h_peticopain-margin, h_petit1copain] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be peticopain
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for peticopain'
    else:
        all_correct = False
        print 'Incorrect result for peticopain. Winning candidate should be 0 and indices are', highest
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        if wa.maxx >= eta1*5/13:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase nocoda
        wa.incnocoda(eta1)
        #decrease alignr
        if wa.alignrC >= eta1:
            wa.incalignrC(-eta1)
        #decrease t1
        if wa.t1 > eta2:
            wa.inct1(-eta2*3/4) #Because we have 4 examples that decrease t1 as opposed to 3 for the others
            print 'decreased t1'


    #########################################################################################################
    #Begin h-words
    eta0 = eta0/2
    eta1 = eta1/2
    eta2 = eta2/2

    print 'Checking deux heros'
    #Since hero never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    #Assume that there is no initial schwa
    h_deuhero = wa.alignlV - wa.onset
    h_deuz1hero = -wa.dep*(1-wa.z1) + wa.maxx*wa.z1 #Doesn't satisfy AlignlV
    h_deuz1_hero = -wa.dep*(1-wa.z1) + wa.maxx*wa.z1  +wa.alignrC + wa.alignlV -wa.nocoda - wa.onset
    print 'harmonies for deuhero'
    print h_deuhero, h_deuz1hero, h_deuz1_hero 
    harmonies = [h_deuhero-margin, h_deuz1hero, h_deuz1_hero]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be deuhero
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for deuhero'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.z1 > eta2:
            wa.incz1(-eta2)
            print 'decreased z1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    print
    print


    print 'Checking petit heros.'
    #Since hero never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    #Assume some kind of schwa activation on the left of heros
    h_petihero = wa.alignlV -wa.onset
    h_petit1hero = -wa.dep*(1-wa.t1) + wa.maxx*wa.t1 
    h_petit1_hero = -wa.dep*(1-wa.t1) + wa.maxx*wa.t1 +wa.alignrC -wa.nocoda - wa.onset + wa.alignlV 
    print 'harmonies for petihero'
    print h_petihero, h_petit1hero, h_petit1_hero 
    harmonies = [h_petihero-margin, h_petit1hero, h_petit1_hero]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be petihero
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for petihero'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.t1 > eta2:
            wa.inct1(-eta2*3/4)#Because we have 4 examples that decrease t1 as opposed to 3 for the others
            print 'decreased t1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    print
    print
    print 'Checking bon heros.'
    #Since hero never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    h_bohero = wa.alignlV - wa.onset
    h_bon1hero = -wa.dep*(1-wa.n1) + wa.maxx*wa.n1 
    h_bon1_hero = -wa.dep*(1-wa.n1) + wa.maxx*wa.n1 +wa.alignlV  +wa.alignrC -wa.nocoda - wa.onset
    print 'harmonies for bohero'
    print h_bohero, h_bon1hero, h_bon1_hero 
    harmonies = [h_bohero-margin, h_bon1hero, h_bon1_hero]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be bohero
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for bohero'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.n1 > eta2:
            wa.incn1(-eta2)
            print 'decreased n1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    print
    print


    #Start hasard with three liaison consonants
    ######################################################################

    print 'Checking petit hasard.'
    #Since hasard never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    h_petiazar =  -wa.onset + wa.alignlV
    h_petit1azar = -wa.dep*(1-wa.t1) + wa.maxx*wa.t1 
    h_petit1_azar = -wa.dep*(1-wa.t1) + wa.maxx*wa.t1 +wa.alignrC -wa.nocoda - wa.onset + wa.alignlV
    print 'harmonies for petiazar'
    print h_petiazar-margin, h_petit1azar, h_petit1_azar 
    harmonies = [h_petiazar, h_petit1azar, h_petit1_azar]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be petiazar
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for petiazar'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.t1 > eta2:
            wa.inct1(-eta2*3/4)#Because we have 4 examples that decrease t1 as opposed to 3 for the others
            print 'decreased t1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    print
    print


    print 'Checking deux hasards.'
    #Since hasard never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    h_deuazar =  -wa.onset + wa.alignlV
    h_deuz1azar = -wa.dep*(1-wa.z1) + wa.maxx*wa.z1 
    h_deuz1_azar = -wa.dep*(1-wa.z1) + wa.maxx*wa.z1 +wa.alignrC -wa.nocoda - wa.onset + wa.alignlV
    print 'harmonies for deuazar'
    print h_deuazar, h_deuz1azar, h_deuz1_azar 
    harmonies = [h_deuazar-margin, h_deuz1azar, h_deuz1_azar]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be deuazar
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for deuazar'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.z1 > eta2:
            wa.incz1(-eta2*3/4)#Because we have 4 examples that decrease z1 as opposed to 3 for the others
            print 'decreased z1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    print
    print


    print 'Checking bon hasard.'
    #Since hasard never occurs with initial t z or n we assume that the learner will not posit these initial consonants in the UR.
    h_boazar =  -wa.onset + wa.alignlV
    h_bon1azar = -wa.dep*(1-wa.n1) + wa.maxx*wa.n1 
    h_bon1_azar = -wa.dep*(1-wa.n1) + wa.maxx*wa.n1 +wa.alignrC -wa.nocoda - wa.onset + wa.alignlV
    print 'harmonies for boazar'
    print h_boazar, h_bon1azar, h_bon1_azar 
    harmonies = [h_boazar-margin, h_bon1azar, h_bon1_azar]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be boazar
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for boazar'
    else:
        all_correct = False
        if 2 not in highest:
            #decrease onset
            if wa.onset > eta1+0.01:
                wa.inconset(-eta1)
        if 1 not in highest:
            #increase alignl
            wa.incalignlV(eta1)
            print 'increasing alignlV.'
        #increase dep
        wa.incdep(eta1*5/13)
        #decrease max
        wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        if wa.n1 > eta2:
            wa.incn1(-eta2*3/4)#Because we have 4 examples that decrease n1 as opposed to 3 for the others
            print 'decreased n1'
        if 2 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)

    #End hasard with three liaison consonants
    #########################################################################

    print 'Checking joli heros.'
    h_jolihero = -wa.onset +wa.alignlV  + wa.alignrV +wa.maxx #This max is for the i segment which we have been ignoring here
    h_jolero =  0
    h_jolitero = wa.alignrV + wa.maxx -wa.dep #There is no t in the input so we have the full value of dep
    h_jolit_ero = wa.maxx -wa.dep +wa.alignlC  -wa.nocoda - wa.onset
    print 'Harmonies for jolihero'
    print h_jolihero, h_jolero, h_jolitero, h_jolit_ero
    harmonies = [h_jolihero - margin, h_jolero, h_jolitero, h_jolit_ero]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be jolihero
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for jolihero'
    else:
        all_correct = False
        print 'Incorrect result for jolihero'
        #Decrease onset
        if wa.onset > eta1:
            wa.inconset(-eta1)
        if 2 not in highest:
            wa.incalignrV(eta1)
        wa.incalignlV
        if 1 in highest:
            wa.incmax(eta1)

    print
    print


    print 'Checking le hasard.'
    #le hasard
    h_le1hasard =  wa.alignrV -wa.dep*(1-wa.ele) + wa.maxx*wa.ele - wa.onset + wa.alignlV
    h_lasard = 0
    print 'Harmonies for le hasard'
    print h_le1hasard, h_lasard
    print 'eha', wa.eha
    harmonies =[ h_le1hasard-margin, h_lasard]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We do not want the last member to be highest.
    if 1 not in highest:
        #We have the correct result
        print 'Correct result for le hasard'
    else:
        all_correct = False
        print 'Incorrect result for le hasard'
        #Decrease onset
        if wa.onset >= eta1:
            wa.inconset(-eta1)
        #decrease dep and increase max
        if wa.dep >= eta1*5/13:
            wa.incdep(-eta1*5/13)
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        wa.incalignlV(eta1)
        wa.incalignrV(eta1)


    print
    print


    print 'Checking le heros.'
    #le hero
    h_le1hero =  wa.alignrV -wa.dep*(1-wa.ele) + wa.maxx*wa.ele - wa.onset + wa.alignlV
    h_lero = 0
    print 'Harmonies for le hero'
    print h_le1hero, h_lero
    print 'ehe', wa.ehe
    harmonies =[ h_le1hero-margin,  h_lero]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We do not want the last member to be highest.
    if 1 not in highest:
        #We have the correct result
        print 'Correct result for le hero'
    else:
        all_correct = False
        print 'Incorrect result for le hero'
        #Decrease onset
        if wa.onset >= eta1:
            wa.inconset(-eta1)
        #decrease dep and increase max
        if wa.dep >= eta1*5/13:
            wa.incdep(-eta1*5/13)
        wa.incmax(eta1*5/13)
        print 'decreasing dep and increasing max'
        #increase alignr and alignl
        wa.incalignrV(eta1)
        wa.incalignlV(eta1)



    print
    print
    #quel hasard
    print 'Now checking what happens with quel hasard if we can get enchainement. Assume full activation on the final C of quel'
    h_quel_asar = wa.alignrC  -wa.nocoda - wa.onset + wa.alignlV
    h_que_lasar = 0
    harmonies = [h_quel_asar,  h_que_lasar] 
    print 'quel_asar (no ench)', h_quel_asar,  'que-lasar (ench)',h_que_lasar 
    if h_quel_asar > h_que_lasar:
        print 'Correct result for quel hasard'
    else:
        all_correct = False
        #In this situation, enchainement is too much favoured so increase inner boundary
        print 'Enchainement is too favoured.'
        #increase alignr
        wa.incalignrC(eta1)
        #Decrease nocoda and onset
        if wa.nocoda >= eta1:
            wa.incnocoda(-eta1)
        if wa.onset >= eta1:
            wa.inconset(-eta1)
        #increase alignl
        wa.incalignlV(eta1)
        print 'increasing alignl.'

    print
    print


    print
    print


    print 'Now checking what happens with quel heros if we get enchainement which we dont want. Assume full activation on the final C of quel'
    h_quel_heros = wa.alignrC  -wa.nocoda - wa.onset + wa.alignlV
    h_que_leros = 0
    harmonies = [h_quel_heros,  h_que_leros] 
    print 'no enchainement', h_quel_heros, 'enchainement',  h_que_leros 

    if h_que_leros + margin  > h_quel_heros:
        all_correct = False
        print 'Wrong result for quel heros with values', round(h_quel_heros,2), 'for non echainement and', round(h_que_leros,2), 'for echainement'
        wa.incalignrC(eta1)
        #Decrease nocoda and onset
        if wa.nocoda >= eta1:
            wa.incnocoda(-eta1)
        if wa.onset >= eta1:
            wa.inconset(-eta1)
        #increase alignl
        wa.incalignlV(eta1)
        print 'increasing alignl.'
    else:
        print 'Correct result for quel heros with no enchainement.'

    print
    print

    print
    print
    #End of h words
    ####################################################################################

    eta0 = eta0*2
    eta1 = eta1*2
    eta2 = eta2*2

    h_petichat = wa.alignlC
    h_petitchat = -wa.dep*(1-wa.t1) + wa.maxx*(wa.t1) - wa.nocoda + wa.alignlC + wa.alignrC
    print 'Harmonies for petichat'
    print h_petichat, h_petitchat
    harmonies =[h_petichat - margin, h_petitchat]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be petichat
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for petichat'
    else:
        all_correct = False
        print 'Incorrect result for petit chat'
        if wa.alignrC > eta1:
            wa.incalignrC(-eta1)
        wa.incdep(eta1)
        if wa.maxx > eta1:
            wa.incmax(-eta1)
        wa.incnocoda(eta1)

    print
    print


    h_deuchat = wa.alignlC
    h_deuzchat = -wa.dep*(1-wa.z1) + wa.maxx*(wa.z1) - wa.nocoda + wa.alignlC + wa.alignrC
    print 'Harmonies for deuchat'
    print h_deuchat, h_deuzchat
    harmonies =[h_deuchat - margin, h_deuzchat]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be deuchat
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for deuchat'
    else:
        all_correct = False
        print 'Incorrect reslt for deux chats'
        if wa.alignrC > eta1:
            wa.incalignrC(-eta1)
        wa.incdep(eta1)
        if wa.maxx > eta1:
            wa.incmax(-eta1)
        wa.incnocoda(eta1)

    print
    print

    h_bochat = wa.alignlC
    h_bonchat = -wa.dep*(1-wa.n1) + wa.maxx*(wa.n1) - wa.nocoda + wa.alignlC + wa.alignrC
    print 'Harmonies for bochat'
    print h_bochat, h_bonchat
    harmonies =[h_bochat - margin, h_bonchat]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be bochat
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for bochat'
    else:
        all_correct = False
        print 'Incorrect resut for bonchat'
        if wa.alignrC > eta1:
            wa.incalignrC(-eta1)
        wa.incdep(eta1)
        if wa.maxx > eta1:
            wa.incmax(-eta1)
        wa.incnocoda(eta1)

    print
    print


    print 'Checking l\'ami.'
    #l'ami
    h_leami = wa.alignrV - wa.dep*(1-wa.ele) + wa.maxx*wa.ele - wa.onset 
    h_lami = 0
    if mult_allo:
        (allomorph, allomorph_name) = wa.chooseL1()
        if allomorph_name == 'n2':
            h_lenami =  wa.alignrV -wa.dep*(1-wa.ele) -wa.dep*(1-wa.n2) + wa.maxx*(wa.ele+wa.n2) + wa.alignlC
            h_len_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.n2) + wa.maxx*(wa.ele+wa.n2) -wa.onset -wa.nocoda 
            h_letami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_let_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
            h_lezami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_lez_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
        if allomorph_name == 't2':
            h_letami =  wa.alignrV -wa.dep*(1-wa.ele) -wa.dep*(1-wa.t2) + wa.maxx*(wa.ele+wa.t2) + wa.alignlC
            h_let_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.t2) + wa.maxx*(wa.ele+wa.t2) -wa.onset -wa.nocoda 
            h_lenami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_len_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
            h_lezami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_lez_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
        if allomorph_name == 'z2':
            h_lezami =  wa.alignrV -wa.dep*(1-wa.ele) -wa.dep*(1-wa.z2) + wa.maxx*(wa.ele+wa.z2) + wa.alignlC
            h_lez_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.z2) + wa.maxx*(wa.ele+wa.z2) -wa.onset -wa.nocoda 
            h_lenami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_len_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
            h_letami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_let_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
        if allomorph_name == 'null':
            h_lezami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) 
            h_lez_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda 
            h_lenami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_len_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
            h_letami =  wa.alignrV -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele)
            h_let_ami =  -wa.dep*(2-wa.ele) + wa.maxx*(wa.ele) -wa.onset -wa.nocoda
    else:
            h_letami =  wa.alignrV -wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[0]) + wa.maxx*(wa.ele+wa.combined[0]) + wa.alignlC
            h_let_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[0]) + wa.maxx*(wa.ele+wa.combined[0]) -wa.onset -wa.nocoda 
            h_lezami =  wa.alignrV -wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[1]) + wa.maxx*(wa.ele+wa.combined[1]) + wa.alignlC
            h_lez_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[1]) + wa.maxx*(wa.ele+wa.combined[1]) -wa.onset -wa.nocoda 
            h_lenami =  wa.alignrV-wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[2]) + wa.maxx*(wa.ele+wa.combined[2]) + wa.alignlC
            h_len_ami =  -wa.dep*(1-wa.ele) -wa.dep*(1-wa.combined[2]) + wa.maxx*(wa.ele+wa.combined[2]) -wa.onset -wa.nocoda 
    print 'Harmonies for l\'ami'
    print h_leami, h_lami, h_lenami, h_len_ami, h_letami, h_let_ami, h_lezami, h_lez_ami
    harmonies = [ h_leami, h_lami - margin, h_lenami, h_len_ami, h_letami, h_let_ami, h_lezami, h_lez_ami]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be lami
    if 1 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for lami'
    else:
        all_correct = False
        print 'Incorrect resut for lami'
        #Decrease epsilon for le
        if wa.ele >= eta2:
            wa.incele(-eta2)
        #Decrease alignr
        if wa.alignrV >= eta1:
            wa.incalignrV(-eta1);
        #Increase dep
        wa.incdep(eta1*5/13)
        #Decrease max
        if wa.maxx >= eta1*5/13:
            wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'
        #increase onset
        wa.inconset(eta1)
        if 2 in highest or 4 in highest or 6 in highest:
            #decrease alignl
            if wa.alignlC >= eta1:
                wa.incalignlC(-eta1)
                print 'Decreasing alignl.'

    print
    print


    print
    print
    print 'Checking quel ami.'
    #quel ami
    h_quelami = 0
    h_quel_ami = -wa.onset + wa.alignrC -wa.nocoda
    #Assume that t1 and z1 and n1 will not occur with quel
    #THe candidates with VC.V syllabilification can't possibly win but are put in for completeness. No learning is incldued for these as false winners since they can't be.
    if mult_allo:
        (allomorph, allomorph_name) = wa.chooseL1()
        if allomorph_name == 't2':
            print 'Allomorph is t2.'
            h_quelt2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrC -wa.nocoda
            h_quelz2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_queln2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_quelt2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda -wa.complex
            h_quelz2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
            h_queln2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
        elif allomorph_name == 'z2':
            print 'Allomorph is z2.'
            h_quelz2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrC-wa.nocoda
            h_quelt2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_queln2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_quelz2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda-wa.complex
            h_quelt2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
            h_queln2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
        elif allomorph_name == 'n2':
            print 'Allomorph is n2.'
            h_queln2ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrC-wa.nocoda
            h_quelt2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_quelz2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_queln2_ami = -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) -wa.onset - wa.nocoda-wa.complex
            h_quelt2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
            h_quelz2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
        else:
            print 'Allomorph is null.'
            h_quelt2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_quelz2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_queln2ami = -wa.dep + wa.alignrC-wa.nocoda
            h_quelt2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
            h_quelz2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
            h_queln2_ami = -wa.dep -wa.onset - wa.nocoda-wa.complex
    else:
        h_quelt2ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC + wa.alignrC-wa.nocoda
        h_quelz2ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC + wa.alignrC-wa.nocoda
        h_queln2ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC + wa.alignrC-wa.nocoda
        h_quelt2_ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) -wa.onset - wa.nocoda-wa.complex
        h_quelz2_ami = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) -wa.onset - wa.nocoda-wa.complex
        h_queln2_ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) -wa.onset - wa.nocoda-wa.complex


    print h_quelami, h_quelt2ami, h_quelz2ami, h_queln2ami, h_quelt2_ami, h_quelz2_ami, h_queln2_ami, h_quel_ami
    harmonies = [h_quelami-margin, h_quelt2ami, h_quelz2ami, h_queln2ami, h_quelt2_ami, h_quelz2_ami, h_queln2_ami, h_quel_ami]
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible member of highest to be quelami
    if 0 in highest and len(highest) == 1:
        #We have the correct result
        print 'Correct result for quelami'
    else:
        all_correct = False
        print 'Wrong result for quelami. Winning index should be 0 and are:', highest
        wa.incnocoda(eta1)
        if 1 in highest or 2 in highest or 3 in highest:
            #decrease alignlC
            if wa.alignlC > eta1:
                wa.incalignlC(-eta1)
                print 'decreasing alignlC.'
            #decrease alignrC
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)
        if 7 in highest:
            #increase onset.
            wa.inconset(eta1)
        if list(set(highest) & set([1,2,3,4,5,6])) != []:
            #increase dep
            wa.incdep(eta1*5/13)
            #decrease max
            if wa.maxx > eta1*5/13:
                wa.incmax(-eta1*5/13)
        print 'increasing dep and decreasing max'

        if 1 in highest:
            #decrease t2
            if mult_allo:
                if wa.t2 > eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2 by', eta2/2
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased combined[0] by', eta2/2, wa.combined[0]

        if 2 in highest:
            #decrease z2
            if mult_allo:
                if wa.z2 > eta2/2:
                    wa.incz2(-eta2/2)
                    print 'decreased z2'
            else:
                if wa.combined[1] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,1)
                    print 'decreased comb1', wa.combined[1]

        if 3 in highest:
            #decrease n2
            if mult_allo:
                if wa.n2 > eta2/2:
                    wa.incn2(-eta2/2)
                    print 'decreased n2'
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]

    
    #cafes extras


    print
    print
    print 'Checking cafes extras'
    h_cafeextra = -wa.onset
    #t1 may start out as zero.
    h_cafez1extra = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs)
    h_cafez1_extra = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that z1 and t1 will not occur with cafez
    h_cafenextra = -wa.dep #Not right anchored as long as z1 is there and also not left anchored
    h_cafen_extra = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as z1 is there 
    h_cafetextra = -wa.dep #Not right anchored as long as z1 is there and also not left anchored
    h_cafet_extra = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as z1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'z2':
            #Input is cafe(z1) + z2extra
            print 'Allomorph is z2.'
            h_cafez2extra =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_cafez2_extra =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_cafez12extra = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignlC - wa.unif*wa.mu
            h_cafez12extraz = -wa.dep*(2.0 - min((wa.zs*2 + allomorph),1)) + wa.maxx*min(1,(wa.zs*2 + allomorph)) + wa.alignlC - wa.unif*wa.mu - wa.nocoda
            h_cafez12_extra = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
            h_cafez1_t2extra = -100.0
            h_cafez1_n2extra = -100.0
            h_cafen2extra = -100.0
            h_cafen2_extra = -100.0
            h_cafet2extra = -100.0
            h_cafet2_extra = -100.0
        elif allomorph_name == 't2':
            #input is cafe(z1) + t2extra
            print 'Allomorph is t2.'
            h_cafez2extra = -100.0
            h_cafez2_extra = -100.0
            h_cafez12extra = -100.0
            h_cafez12extraz = -100.0
            h_cafez12_extra = -100.0
            h_cafez1_t2extra =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.t2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_cafez1_n2extra = -100.0
            h_cafet2extra =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC 
            h_cafet2_extra =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset
            h_cafen2extra = -100.0
            h_cafen2_extra = -100.0
        elif allomorph_name == 'n2':
            #input is cafe(z1) + n2extra
            print 'Allomorph is n2.'
            h_cafez2extra = -100.0
            h_cafez2_extra = -100.0
            h_cafez12extra = -100.0
            h_cafez12extraz = -100.0
            h_cafez12_extra = -100.0
            h_cafez1_t2extra = -100.0
            h_cafez1_n2extra =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.n2) + wa.alignrC - wa.alignlC - wa.nocoda 
            h_cafet2extra = -100.0
            h_cafet2_extra = -100.0
            h_cafen2extra =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC 
            h_cafen2_extra =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
        else:
            #input is cafez(z1) + extra
            print 'Allomorph is null.'
            h_cafez2extra = -100.0
            h_cafez2_extra = -100.0
            h_cafez12extra = -100.0
            h_cafez12extraz = -100.0
            h_cafez12_extra = -100.0
            h_cafez1_n2extra = -100.0
            h_cafez1_t2extra = -100.0
            h_cafen2extra = -100.0
            h_cafen2_extra = -100.0
            h_cafet2extra = -100.0
            h_cafet2_extra = -100.0
    else:
        print 'Blended t2, z2, n2'
        #input is cafe(z1) + [t2,z2,n2]extra
        h_cafez2extra =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC
        h_cafez2_extra =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_cafez12extra = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignlC - wa.unif*wa.mu
        h_cafez12_extra = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset
        h_cafez1_n2extra = -100.0
        h_cafez1_t2extra = -100.0
        h_cafen2extra = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_cafen2_extra = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset
        h_cafet2extra = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_cafet2_extra = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_cafez12extraz = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1.0,(wa.zs + wa.combined[1])) + -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1.0,wa.zs) + wa.alignlC - wa.unif*wa.mu -wa.nocoda

    print 'harmonies for cafez extras'
    #print h_cafez12extraz 
    print h_cafeextra, h_cafez1extra, h_cafez1_extra, h_cafenextra, h_cafen_extra, h_cafetextra, h_cafet_extra, h_cafez2extra, h_cafez2_extra, h_cafez12extra, h_cafez12_extra, h_cafez1_n2extra, h_cafez1_t2extra, h_cafen2extra, h_cafet2extra, h_cafen2_extra, h_cafet2_extra, h_cafez12extraz 

#                0                1                      2            3               4

    harmonies =[ h_cafeextra, h_cafez1extra-margin*2, h_cafez1_extra, h_cafenextra, h_cafen_extra, h_cafetextra, h_cafet_extra, h_cafez2extra-margin*2, h_cafez2_extra, h_cafez12extra-margin*2, h_cafez12_extra, h_cafez1_n2extra, h_cafez1_t2extra, h_cafen2extra, h_cafet2extra, h_cafen2_extra, h_cafet2_extra, h_cafez12extraz ] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible members of highest to be cafez1extra and cafez12extra and cafez2extra and cafez1_extra and cafez12_extra and cafez2_extra
    bad = [0,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18]
    if list(set(highest) & set(bad)) == []:
    #if any of 1,7,9 in highest:
        print 'For cafezextra, highest exclusively includes at least one correct candidate and should be any of 1,7,9', highest
    else:
        all_correct = False
        print 'Wrong result for cafezextra. Winning candidate indices should be exclusive from the set 1,7,90 and are are', highest
        #We want cafez1extra or cafez12extra to win. We can do the following.
        #increase z1
        if wa.zs <= 1-eta2:
            wa.inczs(eta2)
        #increase z2
        if mult_allo:
            wa.incz2(eta2/2)
            print 'increased z2'
        else:
            wa.inccomb(eta2/2,1)
            print 'increased comb1'
        #decrease uniformity
        if wa.unif*wa.mu > eta0:
            wa.incunif(-eta0)

        #decrease dep as long as it doesn't go below zero
        if 0 in highest: 
            #decrease dep
            if wa.dep > eta1*5/13:
                wa.incdep(-eta1*5/13)
            #increase onset
            wa.inconset(eta1*5/13)
            #increase maxx
            wa.incmax(eta1)
            print 'decreasing dep and increasing max'
        if list(set(highest) & set([11,13,15])) != []:
            #decrease n2
            if mult_allo:
                if wa.n2 >= eta2/2:
                    wa.incn2(-eta2/2)
                    print 'decreased n2'
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]
        if list(set(highest) & set([12,14,16])) != []:
            #decrease t2
            if mult_allo:
                if wa.t2 >= eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2'
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased comb0', wa.combined[0]

        if list(set(highest) & set([2,4,6,8,10,11,12,15,16,17,18])) != []:
            wa.incnocoda(eta1)
            wa.inconset(eta1)

    print



    print
    print
    print 'Checking chers amis'
    #Assume the r to start with full activation.
    h_cherami = 0
    h_cher_ami = wa.nocoda -wa.onset +wa.alignrC
    h_cherz1ami = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) -wa.nocoda +wa.alignrC
    h_chez1ami = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) -wa.maxx
    h_cherz1_ami = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) -wa.onset -wa.nocoda + wa.alignrC -wa.complex 
    #*************************Assume that z1 and t1 will not occur with chers
    h_chernami = -wa.dep - wa.nocoda  +wa.alignrC#Not right anchored as long as z1 is there and also not left anchored
    h_chern_ami = -wa.dep  - wa.nocoda - wa.onset -wa.complex #Not right anchored as long as z1 is there 
    h_chertami = -wa.dep -wa.nocoda  +wa.alignrC#Not right anchored as long as z1 is there and also not left anchored
    h_chert_ami = -wa.dep  - wa.nocoda - wa.onset -wa.complex #Not right anchored as long as z1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'z2':
            #Input is cher(z1) + z2ami
            print 'Allomorph is z2.'
            h_cherz2ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC -wa.nocoda  +wa.alignrC
            h_cherz2_ami =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  -wa.nocoda - wa.onset -wa.complex
            h_cherz12ami = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignlC - wa.unif*wa.mu-wa.nocoda +wa.alignrC 
            h_chez12ami = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignlC - wa.unif*wa.mu -wa.maxx
            h_cherz12amiz = -wa.dep*(2.0 - min((wa.zs*2 + allomorph),1)) + wa.maxx*min(1,(wa.zs*2 + allomorph)) + wa.alignlC - wa.unif*wa.mu - wa.nocoda*2 +wa.alignrC
            h_cherz12_ami = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset-wa.complex
            h_cherz1_t2ami = -100.0
            h_cherz1_n2ami = -100.0
            h_chern2ami = -100.0
            h_chern2_ami = -100.0
            h_chert2ami = -100.0
            h_chert2_ami = -100.0
        elif allomorph_name == 't2':
            #input is cher(z1) + t2ami
            print 'Allomorph is t2.'
            h_cherz2ami = -100.0
            h_cherz2_ami = -100.0
            h_cherz12ami = -100.0
            h_chez12ami = -100.0
            h_cherz12amiz = -100.0
            h_cherz12_ami = -100.0
            h_cherz1_t2ami =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.t2) + wa.alignrC + wa.alignlC - wa.nocoda -wa.complex
            h_cherz1_n2ami = -100.0
            h_chert2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC -wa.nocoda  +wa.alignrC
            h_chert2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset -wa.complex
            h_chern2ami = -100.0
            h_chern2_ami = -100.0
        elif allomorph_name == 'n2':
            #input is cher(z1) + n2ami
            print 'Allomorph is n2.'
            h_cherz2ami = -100.0
            h_cherz2_ami = -100.0
            h_cherz12ami = -100.0
            h_chez12ami = -100.0
            h_cherz12amiz = -100.0
            h_cherz12_ami = -100.0
            h_cherz1_t2ami = -100.0
            h_cherz1_n2ami =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.n2) + wa.alignrC - wa.alignlC - wa.nocoda -wa.complex 
            h_chert2ami = -100.0
            h_chert2_ami = -100.0
            h_chern2ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC  -wa.nocoda +wa.alignrC
            h_chern2_ami =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset -wa.complex
        else:
            #input is cherz(z1) + ami
            print 'Allomorph is null.'
            h_cherz2ami = -100.0
            h_cherz2_ami = -100.0
            h_cherz12ami = -100.0
            h_chez12ami = -100.0
            h_cherz12amiz = -100.0
            h_cherz12_ami = -100.0
            h_cherz1_n2ami = -100.0
            h_cherz1_t2ami = -100.0
            h_chern2ami = -100.0
            h_chern2_ami = -100.0
            h_chert2ami = -100.0
            h_chert2_ami = -100.0
    else:
        print 'Blended t2, z2, n2'
        #input is cher(z1) + [t2,z2,n2]ami
        h_cherz2ami =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC -wa.nocoda +wa.alignrC
        h_cherz2_ami =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset -wa.complex
        h_cherz12ami = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignlC - wa.unif*wa.mu -wa.nocoda +wa.alignrC
        h_chez12ami = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignlC - wa.unif*wa.mu -wa.maxx
        h_cherz12_ami = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignrC - wa.unif*wa.mu - wa.nocoda - wa.onset -wa.complex
        h_cherz1_n2ami = -100.0
        h_cherz1_t2ami = -100.0
        h_chern2ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC -wa.nocoda  +wa.alignrC#No reward for alignr since we assume a t1 to be in the input
        h_chern2_ami = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset -wa.complex
        h_chert2ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC -wa.nocoda  +wa.alignrC#No reward for alignr since we assume a t1 to be in the input
        h_chert2_ami = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset -wa.complex
        h_cherz12amiz = -wa.dep*(2.0 - min((wa.zs*2 + wa.combined[1]),1)) + wa.maxx*min(2,(wa.zs*2 + wa.combined[1])) + wa.alignlC - wa.unif*wa.mu -wa.nocoda +wa.alignrC

    print 'harmonies for chers amis'
    #print h_cherz12amiz 
    print h_cherami, h_cherz1ami, h_cherz1_ami, h_chernami, h_chern_ami, h_chertami, h_chert_ami, h_cherz2ami, h_cherz2_ami, h_cherz12ami, h_cherz12_ami, h_cherz1_n2ami, h_cherz1_t2ami, h_chern2ami, h_chert2ami, h_chern2_ami, h_chert2_ami, h_cherz12amiz, h_chez12ami 

#                0                1                      2            3               4

    harmonies =[ h_cherami, h_cherz1ami-margin*2, h_cherz1_ami, h_chernami, h_chern_ami, h_chertami, h_chert_ami, h_cherz2ami-margin*2, h_cherz2_ami, h_cherz12ami-margin*2, h_cherz12_ami, h_cherz1_n2ami, h_cherz1_t2ami, h_chern2ami, h_chert2ami, h_chern2_ami, h_chert2_ami, h_cherz12amiz, h_chez12ami ] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the only possible members of highest to be cherz1ami and cherz12ami and cherz2ami and cherz1_ami and cherz12_ami and cherz2_ami
    bad = [0,2,3,4,5,6,8,10,11,12,13,14,15,16,17,18]
    if list(set(highest) & set(bad)) == []:
    #if any of 1,7,9 in highest:
        print 'For cherzami, highest exclusively includes at least one correct candidate and should be any of 1,7,9', highest
    else:
        all_correct = False
        print 'Wrong result for cherzami. Winning candidate indices should be exclusive from the set 1,7,9 and are are', highest
        #We want cherz1ami or cherz12ami to win. We can do the following.
        #increase z1
        if wa.zs < 1-eta2:
            wa.inczs(eta2)
        #increase z2
        if mult_allo:
            wa.incz2(eta2/2)
            print 'increased z2'
        else:
            wa.inccomb(eta2/2,1)
            print 'increased comb1'
        #decrease uniformity
        if wa.unif*wa.mu > eta0:
            wa.incunif(-eta0)

        #decrease dep as long as it doesn't go below zero
        if 0 in highest: 
            #decrease dep
            if wa.dep > eta1*5/13:
                wa.incdep(-eta1*5/13)
            #increase onset
            wa.inconset(eta1*5/13)
            #increase maxx
            wa.incmax(eta1)
        if list(set(highest) & set([11,13,15])) != []:
            #decrease n2
            if mult_allo:
                if wa.n2 >= eta2/2:
                    wa.incn2(-eta2/2)
                    print 'decreased n2'
            else:
                if wa.combined[2] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,2)
                    print 'decreased comb2', wa.combined[2]
        if list(set(highest) & set([12,14,16])) != []:
            #decrease t2
            if mult_allo:
                if wa.t2 >= eta2/2:
                    wa.inct2(-eta2/2)
                    print 'decreased t2'
            else:
                if wa.combined[0] > (eta2/2)+combined_limit:
                    wa.inccomb(-eta2/2,0)
                    print 'decreased comb0', wa.combined[0]

        if list(set(highest) & set([2,4,6,8,10,11,12,15,16,17])) != []:
            wa.incnocoda(eta1)
            wa.inconset(eta1)

    print



    print
    print

    if mult_allo:
        if allomorph_name == 't2':
            t2d = wa.t2
            z2d = 0
        elif allomorph_name == 'z2':
            t2d = 0
            z2d = wa.z2
        else:
            t2d = 0
            z2d = 0
    else:
        t2d = wa.combined[0]
        z2d = wa.combined[1]

    print 'Checking petits amis'
    h_petiz12ami = -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*wa.mu + wa.alignlC
    h_petit12z12ami = -wa.dep*(1.0-min(wa.t1+t2d,1)) + wa.maxx*(min(wa.t1+t2d,1)) -wa.dep*(1.0-min(wa.z1+z2d,1)) + wa.maxx*(min(wa.z1+z2d,1)) - wa.unif*wa.mu*2 + wa.alignlC - wa.nocoda + wa.alignrC -wa.integrity
    h_petit12ami = -wa.dep*(1.0-min((wa.t1+t2d),1)) + wa.maxx*(min(wa.t1+t2d,1)) - wa.unif*wa.mu + wa.alignlC 
    h_petiz12_ami = -wa.dep*(1.0-min((wa.zs+z2d),1)) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*wa.mu + wa.alignrC - wa.nocoda -wa.onset
    h_petit12_ami = -wa.dep*(1.0-min((wa.t1+t2d),1)) + wa.maxx*(min(wa.t1+t2d,1)) - wa.unif*wa.mu - wa.nocoda -wa.onset + wa.alignrC
    h_petiz1ami = -wa.dep*(1.0-min((wa.zs),1)) + wa.maxx*(min(wa.zs,1))
    h_petit1z12ami = -wa.dep*(1.0-min((wa.t1),1)) + wa.maxx*(min(wa.t1,1))  -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*wa.mu + wa.alignlC  -wa.nocoda + wa.alignrC
    h_petit1ami = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1))
    h_petit1_ami = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) - wa.nocoda -wa.onset + wa.alignrC
    h_petit2ami = -wa.dep*(1.0-min((t2d,1))) + wa.maxx*(min(t2d,1)) + wa.alignlC
    h_pami = wa.maxx*(-3)
    h_peti_z12t2ami = -wa.dep*(1.0-min(t2d,1)) + wa.maxx*(min(t2d,1))  -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*wa.mu + wa.alignlC  -wa.complex
    h_petit1z1ami = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1))  -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda + wa.alignrC
    h_petit1t2ami = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1))  -wa.dep*(1.0-min((t2d,1))) + wa.maxx*(min(t2d,1)) -wa.nocoda + wa.alignlC + wa.alignrC
    h_petiz1t2ami = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))  -wa.dep*(1.0-min((t2d,1))) + wa.maxx*(min(t2d,1)) -wa.complex + wa.alignrC
    h_petit1z12t2ami =-wa.dep*(1.0-min((t2d,1))) + wa.maxx*(min(t2d,1))  -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1))  - wa.unif*wa.mu + wa.alignlC  -wa.complex -wa.nocoda + wa.alignrC
    h_petit1z2t2ami =-wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1))  - wa.unif*wa.mu -wa.complex -wa.nocoda + wa.alignrC
    h_petit12z1ami = -wa.dep*(1.0-min(wa.t1+t2d,1)) + wa.maxx*(min(wa.t1+t2d,1)) -wa.dep*(1.0-min(wa.z1,1)) + wa.maxx*(min(wa.z1,1)) - wa.unif*wa.mu - wa.nocoda + wa.alignrC

    harmonies =[h_petiz12ami, h_petit12z12ami, h_petit12ami, h_petiz12_ami, h_petit12_ami, h_petiz1ami, h_petit1z12ami, h_petit1ami, h_petit1_ami, h_petit2ami, h_pami, h_peti_z12t2ami, h_petit1z1ami, h_petit1t2ami, h_petiz1t2ami, h_petit1z12t2ami, h_petit1z2t2ami, h_petit12z1ami ]

    print 'harmonies for petitsamis'
    print
    print h_petiz12ami, h_petit12z12ami, h_petit12ami, h_petiz12_ami, h_petit12_ami, h_petiz1ami, h_petit1z12ami, h_petit1ami, h_petit1_ami, h_petit2ami, h_pami, h_peti_z12t2ami, h_petit1z1ami, h_petit1t2ami, h_petiz1t2ami, h_petit1z12t2ami, h_petit1z2t2ami, h_petit12z1ami 
    print

    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest or 3  in highest or 5 in highest:
        print 'Correct result for petitsami with winner(s)', highest
    else:
        #all_correct = False
        print 'Wrong result for petitsami. ', highest

    print
    print
    
    #petits chats
    h_peticha = wa.alignlC
    h_petiz1a = -wa.maxx -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) 
    h_peti_t1cha = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.complex
    h_petit1_cha = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.nocoda + wa.alignlC + wa.alignrC
    h_petiz1cha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) + wa.alignlC + wa.alignrC -wa.nocoda
    h_peti_z1cha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex
    h_petit1_z1cha = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex -wa.nocoda + wa.alignrC
    h_petit1z1_cha = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex -wa.nocoda + wa.alignlC + wa.alignrC

    h_petichaz = wa.alignlC -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))
    h_petiz1az = -wa.maxx -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) 
    h_peti_t1chaz = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.complex -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))
    h_petit1_chaz = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.nocoda + wa.alignlC -wa.nocoda  -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))  + wa.alignrC  + wa.alignrC
    h_petiz1chaz = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) + wa.alignlC + wa.alignrC -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda 
    h_peti_z1chaz = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))
    h_petit1_z1chaz = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex -wa.nocoda -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))  + wa.alignrC
    h_petit1z1_chaz = -wa.dep*(1.0-min((wa.t1,1))) + wa.maxx*(min(wa.t1,1)) -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) - wa.complex -wa.nocoda + wa.alignlC + wa.alignrC -wa.nocoda  +wa.alignrC -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1))

    harmonies =[h_peticha, h_petiz1a, h_peti_t1cha, h_petit1_cha, h_petiz1cha, h_peti_z1cha, h_petit1_z1cha, h_petit1z1_cha, h_petichaz, h_petiz1az, h_peti_t1chaz, h_petit1_chaz, h_petiz1chaz, h_peti_z1chaz, h_petit1_z1chaz, h_petit1z1_chaz]
    print 'harmonies for petitchat'
    print
    print  h_peticha, h_petiz1a, h_peti_t1cha, h_petit1_cha, h_petiz1cha, h_peti_z1cha, h_petit1_z1cha, h_petit1z1_cha, h_petichaz, h_petiz1az, h_peti_t1chaz, h_petit1_chaz, h_petiz1chaz, h_peti_z1chaz, h_petit1_z1chaz, h_petit1z1_chaz

    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for petits chats with winner(s)', highest
    else:
        #all_correct = False
        print 'Wrong result for petits chats. ', highest

    print
    print

    #chers chats
    h_chercha = -wa.nocoda + wa.alignlC + wa.alignrC
    h_cherchaz = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda + wa.alignlC + wa.alignrC + wa.alignrC
    h_cherzcha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda + wa.alignlC + wa.alignrC
    h_chezcha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) + wa.alignlC + wa.alignrC - wa.maxx
    h_cher_zcha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda -wa.complex + wa.alignrC
    h_che_zcha = -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) -wa.nocoda -wa.complex - wa.maxx

    harmonies = [h_chercha, h_cherchaz, h_cherzcha, h_chezcha, h_cher_zcha, h_che_zcha]
    print 'Harmonies for chers chats', h_chercha, h_cherchaz, h_cherzcha, h_chezcha, h_cher_zcha, h_che_zcha
    
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for chers chats with winner(s)', highest
    else:
        all_correct = False
        print 'Wrong result for chers chats. ', highest
        if wa.zs > eta2:
            wa.inczs(-eta2)
        if 4 in highest or 5 in highest:
            wa.incalignlC(eta1)
            wa.inccomplex(eta1)
        wa.incdep(eta1)
        if wa.maxx > eta1:
            wa.incmax(-eta1)
        if 3 in highest:
            if wa.nocoda > eta1:
                wa.incnocoda(-eta1)
        if 2 in highest or 3 in highest:
            if wa.alignrC > eta1:
                wa.incalignrC(-eta1)
    print
    print

    #amis All plurals when behaving singly will behave like this
    h_ami = wa.alignrV
    h_amiz = wa.alignrC  - wa.nocoda -wa.dep*(1.0-min((wa.zs,1))) + wa.maxx*(min(wa.zs,1)) 
    print 'Harmonies for amis: ami', h_ami, 'amis', h_amiz
    if h_ami-margin > h_amiz:
        print 'Correct result for amis'
    else:
        print 'Wrong result for amis.'
        all_correct = False
        if wa.alignrC > eta1:
            wa.incalignrC(-eta1)
        wa.incnocoda(eta1)
        if wa.zs > eta2*8:
            wa.inczs(-eta2*8)
        wa.incdep(eta1*2)
        if wa.maxx > eta1:
            wa.incmax(-eta1)
        wa.incalignrV(eta1)
    print
    print



    print 'Checking serait age'
    h_seraiage = -wa.onset
    #t1 may start out as zero.
    h_serait1age = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1)
    h_serait1_age = -wa.dep*(1.0 - min(wa.t1,1)) + wa.maxx*min(1,wa.t1) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that z1 and n1 will not occur with serait
    h_seraizage = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_seraiz_age = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    h_serainage = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_serain_age = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 't2':
            #Input is serai(t1) + t2age
            print 'Allomorph is t2.'
            h_serait2age =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_serait2_age =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_serait12age = -wa.dep*(1.0 - min((wa.t1 + allomorph),1)) + wa.maxx*min(1,(wa.t1 + allomorph)) + wa.alignlC - wa.unif
            h_serait12_age = -wa.dep*(1.0 - min((wa.t1 + allomorph),1)) + wa.maxx*min(1,(wa.t1 + allomorph)) + wa.alignrC - wa.unif - wa.nocoda - wa.onset
            h_serait1_z2age = -100.0
            h_serait1_n2age = -100.0
            h_seraiz2age = -100.0
            h_seraiz2_age = -100.0
            h_serain2age = -100.0
            h_serain2_age = -100.0
        elif allomorph_name == 'z2':
            #input is serai(t1) + z2age
            print 'Allomorph is z2.'
            h_serait2age = -100.0
            h_serait2_age = -100.0
            h_serait12age = -100.0
            h_serait12_age = -100.0
            h_serait1_z2age =  -wa.dep*(1.0-min(wa.t1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t1) + wa.maxx*min(1,wa.z2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_serait1_n2age = -100.0
            h_seraiz2age =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  + wa.alignlC 
            h_seraiz2_age =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z2)  - wa.nocoda - wa.onset
            h_serain2age = -100.0
            h_serain2_age = -100.0
            #There should be no colaseced candidates
        elif allomorph_name == 'n2':
            #input is serai(t1) + n2age
            print 'Allomorph is n2.'
            h_serait2age = -100.0
            h_serait2_age = -100.0
            h_serait12age = -100.0
            h_serait12_age = -100.0
            h_serait1_z2age = -100.0
            h_serait1_n2age =  -wa.dep*(1.0-min(wa.t1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t1) + wa.maxx*min(1,wa.n2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_seraiz2age = -100.0
            h_seraiz2_age = -100.0
            h_serain2age =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC 
            h_serain2_age =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
            #There should be no colaseced candidates
        else:
            #input is serait(t1) + age
            print 'Allomorph is null.'
            h_serait2age = -100.0
            h_serait2_age = -100.0
            h_serait12age = -100.0
            h_serait12_age = -100.0
            h_serait1_z2age = -100.0
            h_serait1_n2age = -100.0
            h_seraiz2age = -100.0
            h_seraiz2_age = -100.0
            h_serain2age = -100.0
            h_serain2_age = -100.0
            #There should be no colaseced candidates
    else:
        print 'Blended t2, z2, n2'
        #input is serai(t1) + [t2,z2,n2]age
        h_serait2age =  -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC
        h_serait2_age =  -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_serait12age = -wa.dep*(1.0 - min((wa.t1 + wa.combined[0]),1)) + wa.maxx*min(1,(wa.t1 + wa.combined[0])) + wa.alignlC - wa.unif
        h_serait12_age = -wa.dep*(1.0 - min((wa.t1 + wa.combined[0]),1)) + wa.maxx*min(1,(wa.t1 + wa.combined[0])) + wa.alignrC - wa.unif - wa.nocoda - wa.onset
        h_serait1_z2age = -100.0
        h_serait1_n2age = -100.0
        h_seraiz2age = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_seraiz2_age = -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_serain2age = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC #No reward for alignr since we assume a t1 to be in the input
        h_serain2_age = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset

    print 'harmonies for seraitage'
    print h_seraiage, h_serait1age, h_serait1_age, h_seraizage, h_seraiz_age, h_serainage, h_serain_age, h_serait2age, h_serait2_age, h_serait12age, h_serait12_age, h_serait1_z2age, h_serait1_n2age, h_seraiz2age, h_serain2age, h_seraiz2_age, h_serain2_age, h_serait1_z2age 

#                0            1            2           3            4

    harmonies = [ h_seraiage, h_serait1age, h_serait1_age, h_seraizage, h_seraiz_age, h_serainage, h_serain_age, h_serait2age, h_serait2_age, h_serait12age, h_serait12_age, h_serait1_z2age, h_serait1_n2age, h_seraiz2age, h_serain2age, h_seraiz2_age, h_serain2_age, h_serait1_z2age, h_serait1_z2age] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the non-liaison candidate to win here
    print 'dep:', wa.dep, 't1:', wa.t1, 't2', wa.t2, 'combined0', wa.combined[0], 'max', wa.maxx, 'alignlC', wa.alignlC, 'unif', wa.unif 
    if 0 in highest and len(highest) == 1:
        print 'Correct result for seraitage', highest
    else:
        all_correct = False
        print 'Wrong result for seraitage.' 
        #For now do no learning. 
    print




    print 'Checking tamis enorme'
    mumed = 0.8
    h_tamienorme = -wa.onset
    #t1 may start out as zero.
    h_tamiz1enorme = -wa.dep*(1.0 - min(wa.z1,1)) + wa.maxx*min(1,wa.z1)
    h_tamiz1_enorme = -wa.dep*(1.0 - min(wa.z1,1)) + wa.maxx*min(1,wa.z1) -wa.onset -wa.nocoda + wa.alignrC 
    #*************************Assume that z1 and n1 will not occur with tamit
    h_tamitenorme = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_tamit_enorme = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    h_taminenorme = -wa.dep #Not right anchored as long as t1 is there and also not left anchored
    h_tamin_enorme = -wa.dep  - wa.nocoda - wa.onset #Not right anchored as long as t1 is there 
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'z2':
            #Input is tami(z1) + tzenorme
            print 'Allomorph is z2.'
            h_tamiz2enorme =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC
            h_tamiz2_enorme =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_tamiz12enorme = -wa.dep*(1.0 - min((wa.z1 + allomorph),1)) + wa.maxx*min(1,(wa.z1 + allomorph)) + wa.alignlC - wa.unif*mumed
            h_tamiz12_enorme = -wa.dep*(1.0 - min((wa.z1 + allomorph),1)) + wa.maxx*min(1,(wa.z1 + allomorph)) + wa.alignrC - wa.unif*mumed - wa.nocoda - wa.onset
            h_tamiz1_t2enorme = -100.0
            h_tamiz1_n2enorme = -100.0
            h_tamit2enorme = -100.0
            h_tamit2_enorme = -100.0
            h_tamin2enorme = -100.0
            h_tamin2_enorme = -100.0
        elif allomorph_name == 't2':
            #input is tami(z1) + t2enorme
            print 'Allomorph is t2.'
            h_tamiz2enorme = -100.0
            h_tamiz2_enorme = -100.0
            h_tamiz12enorme = -100.0
            h_tamiz12_enorme = -100.0
            h_tamiz1_t2enorme =  -wa.dep*(1.0-min(wa.z1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z1) + wa.maxx*min(1,wa.t2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_tamiz1_n2enorme = -100.0
            h_tamit2enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC 
            h_tamit2_enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset
            h_tamin2enorme = -100.0
            h_tamin2_enorme = -100.0
            #There should be no colaseced candidates
        elif allomorph_name == 'n2':
            #input is tami(t1) + n2enorme
            print 'Allomorph is n2.'
            h_tamiz2enorme = -100.0
            h_tamiz2_enorme = -100.0
            h_tamiz12enorme = -100.0
            h_tamiz12_enorme = -100.0
            h_tamiz1_t2enorme = -100.0
            h_tamiz1_n2enorme =  -wa.dep*(1.0-min(wa.z1,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.z1) + wa.maxx*min(1,wa.n2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_tamit2enorme = -100.0
            h_tamit2_enorme = -100.0
            h_tamin2enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC 
            h_tamin2_enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
            #There should be no colaseced candidates
        else:
            #input is tamit(t1) + enorme
            print 'Allomorph is null.'
            h_tamiz2enorme = -100.0
            h_tamiz2_enorme = -100.0
            h_tamiz12enorme = -100.0
            h_tamiz12_enorme = -100.0
            h_tamiz1_t2enorme = -100.0
            h_tamiz1_n2enorme = -100.0
            h_tamit2enorme = -100.0
            h_tamit2_enorme = -100.0
            h_tamin2enorme = -100.0
            h_tamin2_enorme = -100.0
            #There should be no colaseced candidates
    else:
        print 'Blended t2, z2, n2'
        #input is tami(t1) + [t2,z2,n2]enorme
        h_tamiz2enorme =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC
        h_tamiz2_enorme =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_tamiz12enorme = -wa.dep*(1.0 - min((wa.z1 + wa.combined[1]),1)) + wa.maxx*min(1,(wa.z1 + wa.combined[1])) + wa.alignlC - wa.unif*mumed
        print 'h_tamiz12enorme'
        print h_tamiz12enorme,  -wa.dep*(1.0 - min((wa.z1 + wa.combined[1]),1)),  wa.maxx*min(1,(wa.z1 + wa.combined[1])),   wa.alignlC, - wa.unif*mumed
        print
        h_tamiz12_enorme = -wa.dep*(1.0 - min((wa.z1 + wa.combined[1]),1)) + wa.maxx*min(1,(wa.z1 + wa.combined[1])) + wa.alignrC - wa.unif*mumed - wa.nocoda - wa.onset
        h_tamiz1_t2enorme = -100.0
        h_tamiz1_n2enorme = -100.0
        h_tamit2enorme = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC #No reward for alignr since we assume a z1 to be in the input
        h_tamit2_enorme = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_tamin2enorme = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC #No reward for alignr since we assume a z1 to be in the input
        h_tamin2_enorme = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset

    print 'harmonies for tamis enorme'
    print h_tamienorme, h_tamiz1enorme, h_tamiz1_enorme, h_tamitenorme, h_tamit_enorme, h_taminenorme, h_tamin_enorme, h_tamiz2enorme, h_tamiz2_enorme, h_tamiz12enorme, h_tamiz12_enorme, h_tamiz1_t2enorme, h_tamiz1_n2enorme, h_tamit2enorme, h_tamin2enorme, h_tamit2_enorme, h_tamin2_enorme, h_tamiz1_t2enorme 

#                0            1            2           3            4

    harmonies = [ h_tamienorme, h_tamiz1enorme, h_tamiz1_enorme, h_tamitenorme, h_tamit_enorme, h_taminenorme, h_tamin_enorme, h_tamiz2enorme, h_tamiz2_enorme, h_tamiz12enorme, h_tamiz12_enorme, h_tamiz1_t2enorme, h_tamiz1_n2enorme, h_tamit2enorme, h_tamin2enorme, h_tamit2_enorme, h_tamin2_enorme, h_tamiz1_t2enorme, h_tamiz1_t2enorme] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the non-liaison candidate to win here
    print 'dep:', wa.dep, 't1:', wa.t1, 't2', wa.t2, 'combined0', wa.combined[0], 'max', wa.maxx, 'alignlC', wa.alignlC, 'unif', wa.unif 
    if 0 in highest and len(highest) == 1:
        print 'Correct result for tamitenorme', highest
    else:
        #all_correct = False
        print 'Wrong result for tamitenorme.', highest 
        #For now do no learning. 
    print


    print 'Checking momies enormes'
    #They all respect AlignRC at the right edge, either for the final m of the stem or the z of the plural, whichever surfaes, so we leave those rewards out below.
    h_momieenorme = -wa.onset + wa.alignrV
    #t1 may start out as zero.
    h_momiezsenorme = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) + wa.alignrV
    h_momiezs_enorme = -wa.dep*(1.0 - min(wa.zs,1)) + wa.maxx*min(1,wa.zs) -wa.onset -wa.nocoda + wa.alignrC 
    h_momietenorme = -wa.dep
    h_momiet_enorme = -wa.dep  - wa.nocoda - wa.onset  + wa.alignrV
    h_momienenorme = -wa.dep  + wa.alignrV
    h_momien_enorme = -wa.dep  - wa.nocoda - wa.onset
    #If we have multiple allomorphs then for t2, z2 and n2 we only get those values if we happened to choose that allomorph. Otherwise, we get the proportions that we have ended up with as blends.
    if mult_allo:
        if allomorph_name == 'z2':
            #Input is momie(zs) + tzenorme
            print 'Allomorph is z2.'
            h_momiez2enorme =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph) + wa.alignlC + wa.alignrV
            h_momiez2_enorme =  -wa.dep*(1.0 - min(allomorph,1)) + wa.maxx*min(1,allomorph)  - wa.nocoda - wa.onset
            h_momiezsz2enorme = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignlC - wa.unif*mumed + wa.alignrV
            h_momiezsz2_enorme = -wa.dep*(1.0 - min((wa.zs + allomorph),1)) + wa.maxx*min(1,(wa.zs + allomorph)) + wa.alignrC - wa.unif*mumed - wa.nocoda - wa.onset
            h_momiezs_t2enorme = -100.0
            h_momiezs_n2enorme = -100.0
            h_momiet2enorme = -100.0
            h_momiet2_enorme = -100.0
            h_momien2enorme = -100.0
            h_momien2_enorme = -100.0
        elif allomorph_name == 't2':
            #input is momie(zs) + t2enorme
            print 'Allomorph is t2.'
            h_momiez2enorme = -100.0
            h_momiez2_enorme = -100.0
            h_momiezsz2enorme = -100.0
            h_momiezsz2_enorme = -100.0
            h_momiezs_t2enorme =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.t2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_momiezs_n2enorme = -100.0
            h_momiet2enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  + wa.alignlC  + wa.alignrV
            h_momiet2_enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.t2)  - wa.nocoda - wa.onset
            h_momien2enorme = -100.0
            h_momien2_enorme = -100.0
            #There should be no colaseced candidates
        elif allomorph_name == 'n2':
            #input is momie(t1) + n2enorme
            print 'Allomorph is n2.'
            h_momiez2enorme = -100.0
            h_momiez2_enorme = -100.0
            h_momiezsz2enorme = -100.0
            h_momiezsz2_enorme = -100.0
            h_momiezs_t2enorme = -100.0
            h_momiezs_n2enorme =  -wa.dep*(1.0-min(wa.zs,1)) -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.zs) + wa.maxx*min(1,wa.n2) + wa.alignrC + wa.alignlC - wa.nocoda 
            h_momiet2enorme = -100.0
            h_momiet2_enorme = -100.0
            h_momien2enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  + wa.alignlC  + wa.alignrV
            h_momien2_enorme =  -wa.dep*(1.0-min(allomorph,1)) + wa.maxx*min(1,wa.n2)  - wa.nocoda - wa.onset
            #There should be no colaseced candidates
        else:
            #input is momiet(t1) + enorme
            print 'Allomorph is null.'
            h_momiez2enorme = -100.0
            h_momiez2_enorme = -100.0
            h_momiezsz2enorme = -100.0
            h_momiezsz2_enorme = -100.0
            h_momiezs_t2enorme = -100.0
            h_momiezs_n2enorme = -100.0
            h_momiet2enorme = -100.0
            h_momiet2_enorme = -100.0
            h_momien2enorme = -100.0
            h_momien2_enorme = -100.0
            #There should be no colaseced candidates
    else:
        print 'Blended t2, z2, n2'
        #input is momie(t1) + [t2,z2,n2]enorme
        h_momiez2enorme =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1]) + wa.alignlC + wa.alignrV
        h_momiez2_enorme =  -wa.dep*(1.0 - min(wa.combined[1],1)) + wa.maxx*min(1,wa.combined[1])  - wa.nocoda - wa.onset
        h_momiezsz2enorme = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignlC - wa.unif*mumed + wa.alignrV
        print 'h_momiezsz2enorme'
        print h_momiezsz2enorme,  -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)),  wa.maxx*min(1,(wa.zs + wa.combined[1])),  wa.alignlC, - wa.unif*mumed + wa.alignrV
        print
        h_momiezsz2_enorme = -wa.dep*(1.0 - min((wa.zs + wa.combined[1]),1)) + wa.maxx*min(1,(wa.zs + wa.combined[1])) + wa.alignrC - wa.unif*mumed - wa.nocoda - wa.onset
        h_momiezs_t2enorme = -100.0
        h_momiezs_n2enorme = -100.0
        h_momiet2enorme = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0]) + wa.alignlC  + wa.alignrV#No reward for alignr since we assume a zs to be in the input
        h_momiet2_enorme = -wa.dep*(1.0 - min(wa.combined[0],1)) + wa.maxx*min(1,wa.combined[0])  - wa.nocoda - wa.onset
        h_momien2enorme = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2]) + wa.alignlC  + wa.alignrV#No reward for alignr since we assume a zs to be in the input
        h_momien2_enorme = -wa.dep*(1.0 - min(wa.combined[2],1)) + wa.maxx*min(1,wa.combined[2])  - wa.nocoda - wa.onset

    print 'harmonies for momiez enormes'
    print h_momieenorme, h_momiezsenorme, h_momiezs_enorme, h_momietenorme, h_momiet_enorme, h_momienenorme, h_momien_enorme, h_momiez2enorme, h_momiez2_enorme, h_momiezsz2enorme, h_momiezsz2_enorme, h_momiezs_t2enorme, h_momiezs_n2enorme, h_momiet2enorme, h_momien2enorme, h_momiet2_enorme, h_momien2_enorme, h_momiezs_t2enorme 

#                     0            1                2                   3              4

    harmonies = [ h_momieenorme, h_momiezsenorme, h_momiezs_enorme, h_momietenorme, h_momiet_enorme, h_momienenorme, h_momien_enorme, h_momiez2enorme, h_momiez2_enorme, h_momiezsz2enorme, h_momiezsz2_enorme, h_momiezs_t2enorme, h_momiezs_n2enorme, h_momiet2enorme, h_momien2enorme, h_momiet2_enorme, h_momien2_enorme, h_momiezs_t2enorme, h_momiezs_t2enorme] 
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    #We want the liaison candidate to win here
    print 'dep:', wa.dep, 't1:', wa.t1, 't2', wa.t2, 'combined0', wa.combined[0], 'max', wa.maxx, 'alignlC', wa.alignlC, 'unif', wa.unif 
    if 1 in highest or 7 in highest or 9 in highest:
        print 'Correct result for momiezenorme', highest
    else:
        all_correct = False
        print 'Wrong result for momiezenorme.', highest 
        wa.inczs(eta1)
    print


    #Now take a proportion of each of t2, z2 and n2 and put it in combined
    #Take 20% of each and put into combined

    divisor = 5

    if wa.t2 > 0:
        wa.inccomb(wa.t2/divisor,0)
        wa.inct2(-wa.t2/divisor)
    if wa.z2 > 0:
        wa.inccomb(wa.z2/divisor,1)
        wa.incz2(-wa.z2/divisor)
    if wa.n2 > 0:
        wa.inccomb(wa.n2/divisor,2)
        wa.incn2(-wa.n2/divisor)
    tot = wa.combined[0] + wa.combined[1] + wa.combined[2] 

    wa.decrmu(i)



    print
    print


    if all_correct == True:
        print 'comb', wa.combined, 't1', wa.t1, 't2', wa.t2, 'z1', wa.zs, 'z2', wa.z2, 'n1', wa.n1, 'n2', wa.n2, 'dep', wa.dep, 'max', wa.maxx, 'onset', wa.onset, 'nocoda', wa.nocoda, 'unif', wa.unif, 'alignlC', wa.alignlC, 'alignrC', wa.alignrC, 'alignlV', wa.alignlV, 'alignrV', wa.alignrV
        break


print 'comb', wa.combined, 't1', wa.t1, 't2', wa.t2, 'z1', wa.zs, 'z2', wa.z2, 'n1', wa.n1, 'n2', wa.n2, 'dep', wa.dep, 'max', wa.maxx, 'onset', wa.onset, 'nocoda', wa.nocoda, 'unif', wa.unif, 'alignlC', wa.alignlC, 'alignrC', wa.alignrC, 'alignlV', wa.alignlV, 'alignrV', wa.alignrV

print

print 'dep', round(wa.initialdep,2),round(wa.dep,2)
print 'max', round(wa.initialmaxx,2),round(wa.maxx,2)
print 't2', round(wa.initialt2,2),round(wa.t2,2)
print 'z2', round(wa.initialz2,2),round(wa.z2,2)
print 'n2', round(wa.initialn2,2),round(wa.n2,2)
print 'combined[0]', round(wa.initialcombined[0],2),round(wa.combined[0],2)
print 'combined[1]', round(wa.initialcombined[1],2),round(wa.combined[1],2)
print 'combined[2]', round(wa.initialcombined[2],2),round(wa.combined[2],2)
print 't1', round(wa.initialt1,2),round(wa.t1,2)
print 'z1', round(wa.initialz1,2),round(wa.z1,2)
print 'zs', round(wa.initialzs,2),round(wa.zs,2)
print 'beta', round(wa.initialbeta,2),round(wa.beta,2)
print 'beta2', round(wa.initialbeta2,2),round(wa.beta2,2)
print 'omega', round(wa.initialomega,2),round(wa.omega,2)
print 'omega2', round(wa.initialomega2,2),round(wa.omega2,2)
print 'omega3', round(wa.initialomega3,2),round(wa.omega3,2)
print 'phi', round(wa.initialphi,2),round(wa.phi,2)
print 'n1', round(wa.initialn1,2),round(wa.n1,2)
print 't1f', round(wa.initialt1f,2),round(wa.t1f,2)
print 'z1f', round(wa.initialz1f,2),round(wa.z1f,2)
print 'n1f', round(wa.initialn1f,2),round(wa.n1f,2)
print 'ele', round(wa.initialele,2),round(wa.ele,2)
print 'eha', round(wa.initialeha,2),round(wa.eha,2)
print 'eho', round(wa.initialeha,2),round(wa.eho,2)
print 'eham', round(wa.initialeha,2),round(wa.eham,2)
print 'ehe', round(wa.initialehe,2),round(wa.ehe,2)
print 'null', round(wa.initialnull,2),round(wa.null,2)
print 'onset', round(wa.initialonset,2),round(wa.onset,2)
print 'unif', round(wa.initialunif,2),round(wa.unif,2)
print 'mumed', round(mumed,2)
print 'mu', round(wa.mu,2)
print 'alignlC', round(wa.initialalignlC,2),round(wa.alignlC,2)
print 'alignrC', round(wa.initialalignrC,2),round(wa.alignrC,2)
print 'alignlV', round(wa.initialalignlV,2),round(wa.alignlV,2)
print 'alignrV', round(wa.initialalignrV,2),round(wa.alignrV,2)
print 'nocoda', round(wa.initialnocoda,2),round(wa.nocoda,2)
print 'complex', round(wa.initialcomplex,2),round(wa.complex,2)
print


print 'h_peti', round(h_peti,2)
print 'h_petit', round(h_petit, 2)
print

print 'h_deu', round(h_deu,2)
print 'h_deuz', round(h_deuz, 2)
print

print 'h_bo', round(h_bo,2)
print 'h_bon', round(h_bon, 2)
print

print 'h_peti', round(h_peti,2)
print 'h_petite', round(h_petite, 2)
print

print 'h_dou', round(h_dou,2)
print 'h_douze', round(h_douze, 2)
print

print 'h_bo', round(h_bo,2)
print 'h_bonne', round(h_bonne, 2)
print

print 'h_ami', round(h_ami,2)
print 'h_t2ami', round(h_t2ami,2)
print 'h_z2ami', round(h_z2ami,2)
print 'h_n2ami', round(h_n2ami,2)
print

print 'h_petiami', round(h_petiami,2)
print 'h_petit1ami', round(h_petit1ami,2)
print 'h_petit1_ami', round(h_petit1_ami,2)
print 'h_petizami', round(h_petizami,2)
print 'h_petiz_ami', round(h_petiz_ami,2)
print 'h_petinami', round(h_petinami,2)
print 'h_petin_ami', round(h_petin_ami,2)
print 'h_petit2ami', round(h_petit2ami,2)
print 'h_petit2_ami', round(h_petit2_ami,2)
print 'h_petit12ami', round(h_petit12ami,2)
print 'h_petit12_ami', round(h_petit12_ami,2)
print 'h_petit1_z2ami', round(h_petit1_z2ami,2)
print 'h_petit1_n2ami', round(h_petit1_n2ami,2)
print 'h_petiz2ami', round(h_petiz2ami,2)
print 'h_petin2ami', round(h_petin2ami,2)
print 'h_petiz2_ami', round(h_petiz2_ami,2)
print 'h_petin2_ami', round(h_petin2_ami,2)
print 'h_petit1_z2ami ', round(h_petit1_z2ami ,2)
print


print 'h_boami', round(h_boami,2)
print 'h_bon1ami', round(h_bon1ami,2)
print 'h_bon1_ami', round(h_bon1_ami,2)
print 'h_botami', round(h_botami,2)
print 'h_bot_ami', round(h_bot_ami,2)
print 'h_bozami', round(h_bozami,2)
print 'h_boz_ami', round(h_boz_ami,2)
print 'h_bon2ami', round(h_bon2ami,2)
print 'h_bon2_ami', round(h_bon2_ami,2)
print 'h_bon12ami', round(h_bon12ami,2)
print 'h_bon12_ami', round(h_bon12_ami,2)
print 'h_bon1_z2ami', round(h_bon1_z2ami,2)
print 'h_bon1_t2ami', round(h_bon1_t2ami,2)
print 'h_bot2ami', round(h_bot2ami,2)
print 'h_boz2ami', round(h_boz2ami,2)
print 'h_bot2_ami', round(h_bot2_ami,2)
print 'h_boz2_ami', round(h_boz2_ami,2)
print

print 'h_deuami', round(h_deuami,2)
print 'h_deuz1ami', round(h_deuz1ami,2)
print 'h_deuz1_ami', round(h_deuz1_ami,2)
print 'h_deunami', round(h_deunami,2)
print 'h_deun_ami', round(h_deun_ami,2)
print 'h_deutami', round(h_deutami,2)
print 'h_deut_ami', round(h_deut_ami,2)
print 'h_deuz2ami', round(h_deuz2ami,2)
print 'h_deuz2_ami', round(h_deuz2_ami,2)
print 'h_deuz12ami', round(h_deuz12ami,2)
print 'h_deuz12_ami', round(h_deuz12_ami,2)
print 'h_deuz1_n2ami', round(h_deuz1_n2ami,2)
print 'h_deuz1_t2ami', round(h_deuz1_t2ami,2)
print 'h_deun2ami', round(h_deun2ami,2)
print 'h_deut2ami', round(h_deut2ami,2)
print 'h_deun2_ami', round(h_deun2_ami,2)
print 'h_deut2_ami ', round(h_deut2_ami ,2)
print


print 'h_quelami', round(h_quelami,2)
print 'h_quelt2ami', round(h_quelt2ami,2)
print 'h_quelz2ami', round(h_quelz2ami,2)
print 'h_queln2ami', round(h_queln2ami,2)
print 'h_quelt2_ami', round(h_quelt2_ami,2)
print 'h_quelz2_ami', round(h_quelz2_ami,2)
print 'h_queln2_ami', round(h_queln2_ami,2)
print 'h_quel_ami', round(h_quel_ami,2)
print

print 'h_quel_asar', round(h_quel_asar,2)
print 'h_que_lasar', round(h_que_lasar,2)
print


print 'h_quel_heros', round(h_quel_heros,2)
print 'h_que_leros', round(h_que_leros,2)
print

print 'h_joliami', round(h_joliami,2)
print 'h_jolit2ami', round(h_jolit2ami,2)
print 'h_joliz2ami', round(h_joliz2ami,2)
print 'h_jolin2ami', round(h_jolin2ami,2)
print 'h_jolit2_ami', round(h_jolit2_ami,2)
print 'h_joliz2_ami', round(h_joliz2_ami,2)
print 'h_jolin2_ami', round(h_jolin2_ami,2)
print


print 'h_peticopain', round(h_peticopain,2)
print 'h_petit1copain', round(h_petit1copain,2)
print

print 'h_petit1fcopine', round(h_petit1fcopine,2)
print 'h_peticopine', round(h_peticopine,2)
print

print 'h_deucopain', round(h_deucopain,2)
print 'h_deuzcopain', round(h_deuz1copain,2)
print

print 'h_douz1fcopain', round(h_douz1fcopain,2)
print 'h_doucopain', round(h_doucopain,2)
print

print 'h_bocopain', round(h_bocopain,2)
print 'h_bon1copain ', round(h_bon1copain ,2)
print

print 'h_bocopine', round(h_bocopine,2)
print 'h_bon1fcopine ', round(h_bon1fcopine ,2)
print

print 'h_joliami', round(h_joliami,2)
print 'h_jolit2ami', round(h_jolit2ami,2)
print 'h_joliz2ami', round(h_joliz2ami,2)
print 'h_jolin2ami', round(h_jolin2ami,2)
print 'h_jolit2_ami', round(h_jolit2_ami,2)
print 'h_joliz2_ami', round(h_joliz2_ami,2)
print 'h_jolin2_ami', round(h_jolin2_ami,2)
print


print 'h_petihero', round(h_petihero,2)
print 'h_petit1hero', round(h_petit1hero,2)
print ' h_petit1_hero', round( h_petit1_hero,2)
print


print 'h_deuhero', round(h_deuhero,2)
print 'h_deuz1hero', round(h_deuz1hero,2)
print 'h_deuz1_hero', round(h_deuz1_hero,2)
print

print 'h_bohero', round(h_bohero,2)
print 'h_bon1hero', round(h_bon1hero,2)
print 'h_bon1_hero', round(h_bon1_hero,2)
print

print 'h_jolihero', round(h_jolihero,2)
print 'h_jolero', round(h_jolero,2)
print 'h_jolitero', round(h_jolitero,2)
print 'h_jolit_ero', round(h_jolit_ero,2)
print


print 'h_petichat', round(h_petichat,2)
print 'h_petitchat', round(h_petitchat,2)
print

print 'h_deuchat', round(h_deuchat,2)
print 'h_deuzchat', round(h_deuzchat,2)
print

print 'h_bochat', round(h_bochat,2)
print 'h_bonchat', round(h_bonchat,2)
print

print 'h_leami', round(h_leami,2)
print 'h_lami', round(h_lami,2)
print 'h_lenami', round(h_lenami,2)
print 'h_len_ami', round(h_len_ami,2)
print 'h_letami', round(h_letami,2)
print 'h_let_ami', round(h_let_ami,2)
print 'h_lezami', round(h_lezami,2)
print 'h_lez_ami', round(h_lez_ami,2)
print


print 'h_le1hasard', round(h_le1hasard,2)
print 'h_lasard', round(h_lasard,2)
print


print 'h_petiazar', round(h_petiazar,2)
print 'h_petit1azar',round(h_petit1azar,2)
print 'h_petit1_azar',round(h_petit1_azar,2)
print

print 'h_deuazar', round(h_deuazar,2)
print 'h_deuz1azar',round(h_deuz1azar,2)
print 'h_deuz1_azar',round(h_deuz1_azar,2)
print

print 'h_boazar', round(h_boazar,2)
print 'h_bon1azar',round(h_bon1azar,2)
print 'h_bon1_azar',round(h_bon1_azar,2)
print


print 'h_lehero', round(h_le1hero,2)
print 'h_lero', round(h_lero,2)
print


print 'h_cafeextra', round(h_cafeextra,2)
print 'h_cafez1extra', round(h_cafez1extra,2)
print 'h_cafez1_extra', round(h_cafez1_extra,2)
print 'h_cafenextra', round(h_cafenextra,2)
print 'h_cafen_extra', round(h_cafen_extra,2)
print 'h_cafetextra', round(h_cafetextra,2)
print 'h_cafet_extra', round(h_cafet_extra,2)
print 'h_cafez2extra', round(h_cafez2extra,2)
print 'h_cafez2_extra', round(h_cafez2_extra,2)
print 'h_cafez12extra', round(h_cafez12extra,2)
print 'h_cafez12_extra', round(h_cafez12_extra,2)
print 'h_cafez1_n2extra', round(h_cafez1_n2extra,2)
print 'h_cafez1_t2extra', round(h_cafez1_t2extra,2)
print 'h_cafen2extra', round(h_cafen2extra,2)
print 'h_cafet2extra', round(h_cafet2extra,2)
print 'h_cafen2_extra', round(h_cafen2_extra,2)
print 'h_cafet2_extra ', round(h_cafet2_extra ,2)
print 'h_cafez12extraz', round(h_cafez12extraz,2)
print


print 'h_cherami', round(h_cherami,2)
print 'h_cherz1ami', round(h_cherz1ami,3)
print 'h_cherz1_ami', round(h_cherz1_ami,2)
print 'h_chernami', round(h_chernami,2)
print 'h_chern_ami', round(h_chern_ami,2)
print 'h_chertami', round(h_chertami,2)
print 'h_chert_ami', round(h_chert_ami,2)
print 'h_cherz2ami', round(h_cherz2ami,2)
print 'h_cherz2_ami', round(h_cherz2_ami,2)
print 'h_cherz12ami', round(h_cherz12ami,3)
print 'h_chez12ami', round(h_chez12ami,2)
print 'h_cherz12_ami', round(h_cherz12_ami,2)
print 'h_cherz1_n2ami', round(h_cherz1_n2ami,2)
print 'h_cherz1_t2ami', round(h_cherz1_t2ami,2)
print 'h_chern2ami', round(h_chern2ami,2)
print 'h_chert2ami', round(h_chert2ami,2)
print 'h_chern2_ami', round(h_chern2_ami,2)
print 'h_chert2_ami ', round(h_chert2_ami ,2)
print 'h_cherz12amiz', round(h_cherz12amiz,2)
print


print 'h_petiz12ami', round(h_petiz12ami,2)
print 'h_petit12z12ami', round(h_petit12z12ami,2)
print 'h_petit12ami', round(h_petit12ami,2)
print 'h_petiz12_ami', round(h_petiz12_ami,2)
print 'h_petit12_ami', round(h_petit12_ami,2)
print 'h_petiz1ami', round(h_petiz1ami,2)
print 'h_petit1z12ami', round(h_petit1z12ami,2)
print 'h_petit1ami', round(h_petit1ami,2)
print 'h_petit1_ami', round(h_petit1_ami,2)
print 'h_petit2ami', round(h_petit2ami,2)
print 'h_pami', round(h_pami,2)
print 'h_peti_z12t2ami', round(h_peti_z12t2ami,2)
print 'h_petit1z1ami', round(h_petit1z1ami,2)
print 'h_petit1t2ami', round(h_petit1t2ami,2)
print 'h_petiz1t2ami', round(h_petiz1t2ami,2)
print 'h_petit1z12t2ami', round(h_petit1z12t2ami,2)
print 'h_petit1z2t2ami', round(h_petit1z2t2ami,2)
print
    
print 'h_peticha', round(h_peticha,2)
print 'h_petiz1a', round(h_petiz1a,2)
print 'h_peti_t1cha', round(h_peti_t1cha,2)
print 'h_petit1_cha', round(h_petit1_cha,2)
print 'h_petiz1cha', round(h_petiz1cha,2)
print 'h_peti_z1cha', round(h_peti_z1cha,2)
print 'h_petit1_z1cha', round(h_petit1_z1cha,2)
print 'h_petit1z1_cha', round(h_petit1z1_cha,2)
print 'h_petichaz', round(h_petichaz,2)
print 'h_petiz1az', round(h_petiz1az,2)
print 'h_peti_t1chaz', round(h_peti_t1chaz,2)
print 'h_petit1_chaz', round(h_petit1_chaz,2)
print 'h_petiz1chaz', round(h_petiz1chaz,2)
print 'h_peti_z1chaz', round(h_peti_z1chaz,2)
print 'h_petit1_z1chaz', round(h_petit1_z1chaz,2)
print 'h_petit1z1_chaz', round(h_petit1z1_chaz,2)


#Keep all the values learned above the same.

for i in range(50):
    all_correct = True
    print 'iteration', i
    print

    #beau copain All candidates repsect alignreither V or C
    h_bokope = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity 
    h_beokope = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu
    h_be_okope = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset + wa.alignrV
    h_belkope = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.nocoda + wa.alignrC
    harmonies = [h_bokope, h_beokope, h_be_okope, h_belkope]
    print  h_bokope, h_beokope, h_be_okope, h_belkope
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest and len(highest) ==1:
        print 'Correct result for beau copain', highest
    else:   
        all_correct = False
        print 'Wrong result for beau copain.', highest 
        if wa.omega < 1-eta1:
            wa.incomega(eta1)
        if wa.contiguity > eta1:
            wa.inccontig(-eta1)
        if wa.beta > eta1:
            wa.incbeta(-eta1)
        if 1 in highest:
            wa.inccomplex(eta1)


    #beaux copains        
    h_bokope = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity 
    h_bozkope = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) + wa.alignrC -wa.contiguity -wa.nocoda
    h_beokope = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu
    h_be_okope = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset + wa.alignrV
    h_belkope = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.nocoda + wa.alignrC
    h_belzkope = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  -wa.nocoda + wa.alignrC -wa.complex
    harmonies = [h_bokope, h_bozkope, h_beokope, h_be_okope, h_belkope, h_belzkope]
    print h_bokope, h_bozkope, h_beokope, h_be_okope, h_belkope, h_belzkope
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest and len(highest) ==1:
        print 'Correct result for beau copain', highest
    else:   
        all_correct = False
        print 'Wrong result for beau copain.', highest 
        if wa.omega < 1-eta1:
            wa.incomega(eta1)
        if wa.contiguity > eta1:
            wa.inccontig(-eta1)
        if wa.beta > eta1:
            wa.incbeta(-eta1)
        if 1 in highest or 2 in highest:
            wa.inccomplex(eta1)

    #belle copine
    h_belfkopin = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega + wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) -wa.nocoda + wa.alignrC
    h_beflkopin = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.nocoda + wa.alignrC
    h_belkopin = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.nocoda + wa.alignrC
    h_bokopin = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity 
    h_bofkopin = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) + wa.alignrV -wa.contiguity
    h_beokopin = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu
    h_be_okopin = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset + wa.alignrV
    harmonies = [h_belfkopin, h_beflkopin, h_belkopin, h_bokopin, h_beokopin, h_be_okopin]
    print h_belfkopin, h_beflkopin, h_belkopin, h_bokopin, h_beokopin, h_be_okopin
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest or 1 in highest or 2 in highest:
        print 'Correct result for belle copine', highest
    else:   
        all_correct = False
        print 'Wrong result for belle copine.', highest 
        if 3 in highest or 4 in highest:
            wa.incbeta(eta1)
            wa.inccontig(eta1)
        if 5 in highest:
            wa.inccomplex(eta1)
        if 4 not in highest:
            wa.incphi(eta1)
            

    #bel ami
    h_belami =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))
    h_bel_ami =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrC -wa.nocoda -wa.onset
    h_boami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity -wa.onset 
    h_beoami = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu -wa.onset
    h_be_oami = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset*2 + wa.alignrV
    h_botami = -wa.dep*(1.0-min(wa.omega,1)) -wa.dep*(1.0-min(wa.combined[0],1)) + wa.maxx*(min(wa.omega,1)) + wa.maxx*(min(wa.combined[0],1)) + wa.alignrV -wa.contiguity +wa.alignlC
    harmonies = [h_belami, h_bel_ami, h_boami, h_beoami, h_be_oami, h_botami]
    print h_belami, h_bel_ami, h_boami, h_beoami, h_be_oami, h_botami
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for bel ami', highest
    else:   
        all_correct = False
        print 'Wrong result for bel ami.', highest 
        if 2 in highest or 5 in highest:
            wa.incbeta(eta1)
            wa.inccontig(eta1)
        if 3 in highest:
            wa.inccomplex(eta1)


    #belle amie
    h_belami =  -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))
    h_bel_ami =  -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrC -wa.nocoda -wa.onset
    h_boami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) + wa.alignrV -wa.contiguity -wa.onset 
    h_beoami = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu -wa.onset
    h_be_oami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset*2 + wa.alignrV
    h_botami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) -wa.dep*(1.0-min(wa.combined[0],1)) + wa.maxx*(min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.combined[0],1)) + wa.alignrV -wa.contiguity +wa.alignlC
    harmonies = [h_belami, h_bel_ami, h_boami, h_beoami, h_be_oami, h_botami]
    print h_belami, h_bel_ami, h_boami, h_beoami, h_be_oami, h_botami
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for belle ami', highest
    else:   
        all_correct = False
        print 'Wrong result for belle ami.', highest 
        if 2 in highest or 5 in highest:
            wa.incbeta(eta1)
            wa.inccontig(eta1)
        if 3 in highest:
            wa.inccomplex(eta1)

    #beau heros
    h_boero =  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity -wa.onset +wa.alignlV
    h_belero =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))
    h_bel_ero =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))+ wa.alignlV + wa.alignrC -wa.nocoda -wa.onset
    h_beoero = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu -wa.onset + wa.alignlV
    h_be_oero = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset*2 + wa.alignrV
    harmonies = [h_bohero, h_belero, h_bel_ero, h_beoero, h_be_oero]
    print  h_bohero, h_belero, h_bel_ero, h_beoero, h_be_oero
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for beau heros', highest
    else:   
        all_correct = False
        print 'Wrong result for beau heros.', highest 
        if wa.beta > eta1:
            wa.incbeta(-eta1)
        if wa.contiguity > eta1:
            wa.inccontig(-eta1)


    #beaux amis
    h_bozpami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.zs+wa.combined[1],1)) + wa.maxx*(min(wa.zs+wa.combined[1],1)) -wa.contiguity  +wa.alignlC + wa.alignrV
    h_belami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))
    h_bel_ami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) +wa.alignrC -wa.nocoda -wa.onset
    h_belz_ami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  + wa.alignrC -wa.nocoda -wa.onset -wa.complex
    h_belzami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  + wa.alignrC -wa.nocoda 
    h_belzp2ami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.zs+wa.combined[1],1)) + wa.maxx*(min(wa.zs+wa.combined[1],1))  + wa.alignrC -wa.nocoda +wa.alignlC
    h_boami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity -wa.onset 
    h_botami = -wa.dep*(1.0-min(wa.omega,1)) -wa.dep*(1.0-min(wa.combined[0],1)) + wa.maxx*(min(wa.omega,1)) + wa.maxx*(min(wa.combined[0],1)) + wa.alignrV -wa.contiguity +wa.alignlC
    h_bozami = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) -wa.contiguity + wa.alignrV +wa.alignlC
    harmonies = [h_bozpami, h_belami, h_bel_ami, h_belz_ami, h_belz_ami, h_belzp2ami, h_boami, h_botami, h_bozami]
    print h_bozpami, h_belami, h_bel_ami, h_belz_ami, h_belzami, h_belzp2ami, h_boami, h_botami, h_bozami
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for beaux amis', highest
    else:   
        all_correct = False
        print 'Wrong result for beaux amis.', highest 



    #belles amies
    h_bozpami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) -wa.dep*(1.0-min(wa.zs+wa.combined[1],1)) + wa.maxx*(min(wa.zs+wa.combined[1],1)) -wa.contiguity  +wa.alignlC + wa.alignrV
    h_belami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))
    h_bel_ami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) +wa.alignrC -wa.nocoda -wa.onset
    h_belz_ami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  + wa.alignrC -wa.nocoda -wa.onset -wa.complex
    h_belzami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  + wa.alignrC -wa.nocoda 
    h_belzp2ami = -wa.dep*(1.0-min(wa.beta+wa.phi,1)) + wa.maxx*(min(wa.beta+wa.phi,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))  -wa.dep*(1.0-min(wa.zs+wa.combined[1],1)) + wa.maxx*(min(wa.zs+wa.combined[1],1))  + wa.alignrC -wa.nocoda +wa.alignlC
    h_boami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) + wa.alignrV -wa.contiguity -wa.onset 
    h_botami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) -wa.dep*(1.0-min(wa.combined[0],1)) + wa.maxx*(min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.combined[0],1)) + wa.alignrV -wa.contiguity +wa.alignlC
    h_bozami = -wa.dep*(1.0-min(wa.omega+wa.phi,1)) + wa.maxx*(min(wa.omega+wa.phi,1)) -wa.dep*(1.0-min(wa.zs+wa.combined[1],1)) + wa.maxx*(min(wa.zs+wa.combined[1],1)) -wa.contiguity + wa.alignrV +wa.alignlC
    harmonies = [h_bozpami, h_belami, h_bel_ami, h_belz_ami, h_belz_ami, h_belzp2ami, h_boami, h_botami, h_bozami]
    print h_bozpami, h_belami, h_bel_ami, h_belz_ami, h_belzami, h_belzp2ami, h_boami, h_botami, h_bozami
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 5 in highest:
        print 'Correct result for belles amies', highest
    else:   
        all_correct = False
        print 'Wrong result for belles amies.', highest 



    #beaux heros
    h_boero =  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrV -wa.contiguity -wa.onset +wa.alignlV
    h_belero =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))
    h_bel_ero =  -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) + wa.alignrC -wa.nocoda -wa.onset +wa.alignlV
    h_beoero = -wa.dep*(1.0-min(wa.beta+wa.omega,1)) + wa.maxx*(min(wa.beta+wa.omega,1)) -wa.complex +wa.alignrV -wa.unif*wa.mu -wa.onset + wa.alignlV
    h_be_oero = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1))   -wa.onset*2 + wa.alignrV

    h_belz_ero = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))  + wa.alignrC -wa.nocoda -wa.onset -wa.complex +wa.alignlV
    h_belzero = -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1))  -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))
    h_bozero =  -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) - wa.contiguity + wa.alignrV
    h_be_ozero = -wa.dep*(1.0-min(wa.beta,1)) + wa.maxx*(min(wa.beta,1)) -wa.dep*(1.0-min(wa.omega,1)) + wa.maxx*(min(wa.omega,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1))   -wa.onset + wa.alignrV
    harmonies = [h_bohero, h_belero, h_bel_ero, h_beoero, h_be_oero, h_belz_ero, h_belzero, h_bozero, h_be_ozero]
    print  h_bohero, h_belero, h_bel_ero, h_beoero, h_be_oero, h_belz_ero, h_belzero, h_bozero, h_be_ozero
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for beaux heros', highest
    else:   
        all_correct = False
        print 'Wrong result for beaux heros.', highest 
        if wa.beta > eta1:
            wa.incbeta(-eta1)
        if wa.contiguity > eta1:
            wa.inccontig(-eta1)


    #petites amies

    print 'Checking petites amies'
    h_petit1fzs2ami = -wa.dep*(1.0-min(wa.t1+wa.phi,1)) + wa.maxx*(min(wa.t1+wa.phi,1)) -wa.dep*(1.0-min(wa.zs+z2d,1)) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*(wa.mu+1) + wa.alignlC - wa.nocoda + wa.alignrC
    h_petit1fzs2amiz = -wa.dep*(1.0-min(wa.t1+wa.phi,1)) + wa.maxx*(min(wa.t1+wa.phi,1)) -wa.dep*(1.0-min(wa.zs+z2d,1)) + wa.maxx*(min(wa.zs+z2d,1)) -wa.dep*(1.0-min(z2d,1)) + wa.maxx*(min(z2d,1)) - wa.unif*(wa.mu+1) + wa.alignlC - wa.nocoda*2 + wa.alignrC*2
    h_petit1fzsami = -wa.dep*(1.0-min(wa.t1+wa.phi,1)) + wa.maxx*(min(wa.t1+wa.phi,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) - wa.unif*wa.mu + wa.alignlC - wa.nocoda + wa.alignrC
    h_petit12zsami = -wa.dep*(1.0-min(wa.t1+t2d,1)) + wa.maxx*(min(wa.t1+t2d,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) - wa.unif + wa.alignlC - wa.nocoda + wa.alignrC - wa.contiguity
    h_petit12zs2ami = -wa.dep*(1.0-min(wa.t1+t2d,1)) + wa.maxx*(min(wa.t1+t2d,1)) -wa.dep*(1.0-min(wa.zs+z2d,1)) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif*2 + wa.alignlC - wa.nocoda + wa.alignrC - wa.contiguity
    h_petit1zs2ami = -wa.dep*(1.0-min(wa.t1,1)) + wa.maxx*(min(wa.t1,1)) -wa.dep*(1.0-min(wa.zs+z2d,1)) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif + wa.alignlC - wa.nocoda + wa.alignrC
    h_petizs2ami = -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif + wa.alignlC
    h_petizs2fami = -wa.dep*(1.0-min((wa.zs+z2d+wa.phi,1))) + wa.maxx*(min(wa.zs+z2d+wa.phi,1)) - wa.unif + wa.alignlC
    h_petizs2_ami = -wa.dep*(1.0-min((wa.zs+z2d,1))) + wa.maxx*(min(wa.zs+z2d,1)) - wa.unif - wa.nocoda -wa.onset + wa.alignrC
    h_petizs2f_ami = -wa.dep*(1.0-min((wa.zs+z2d+wa.phi,1))) + wa.maxx*(min(wa.zs+z2d+wa.phi,1)) - wa.unif + wa.alignrC -wa.nocoda -wa.onset
    h_petiami = -wa.onset


    harmonies = [h_petit1fzs2ami, h_petit1fzs2amiz, h_petit1fzsami, h_petit12zsami, h_petit12zs2ami, h_petit1zs2ami, h_petizs2ami, h_petizs2fami, h_petizs2_ami, h_petizs2f_ami, h_petiami]

    print 'harmonies for petites amies'
    print
    print  h_petit1fzs2ami, h_petit1fzs2amiz, h_petit1fzsami, h_petit12zsami, h_petit12zs2ami, h_petit1zs2ami, h_petizs2ami, h_petizs2fami, h_petizs2_ami, h_petizs2f_ami, h_petiami
    print

    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest or 2 in highest or 3 in highest or 4 in highest or 5 in highest:
        print 'Correct result for petites amies with winner(s)', highest
    else:
        all_correct = False
        print 'Wrong result for petites amies. ', highest
        wa.incphi(eta1)

    print
    print





    #liberal
    h_liberal =  -wa.dep*(1.0-min(wa.beta2,1)) + wa.maxx*(min(wa.beta2,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) +wa.alignrC -wa.nocoda
    h_libero =  -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) +wa.alignrV -wa.contiguity
    harmonies = [h_liberal, h_libero]
    print  h_liberal, h_libero
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for liberal', highest
    else:   
        all_correct = False
        print 'Wrong result for liberal.', highest 
        wa.incbeta2(eta1)
        wa.inccontig(eta1)

    #liberale
    h_liberal =  -wa.dep*(1.0-min(wa.beta2+wa.phi,1)) + wa.maxx*(min(wa.beta2+wa.phi,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) +wa.alignrC -wa.nocoda -wa.unif*wa.mu
    h_libero =  -wa.dep*(1.0-min(wa.omega2+wa.phi,1)) + wa.maxx*(min(wa.omega2+wa.phi,1)) +wa.alignrV -wa.contiguity -wa.unif*wa.mu
    harmonies = [h_liberal, h_libero]
    print  h_liberal, h_libero
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for liberale', highest
    else:   
        all_correct = False
        print 'Wrong result for liberale.', highest 
        wa.incbeta2(eta1)
        wa.inccontig(eta1)


    #liberaux
    h_liberal =  -wa.dep*(1.0-min(wa.beta2,1)) + wa.maxx*(min(wa.beta2,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) +wa.alignrC -wa.nocoda
    h_liberalo =  -wa.dep*(1.0-min(wa.beta2,1)) + wa.maxx*(min(wa.beta2,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) -wa.dep*(1.0-min(wa.omega3,1)) + wa.maxx*(min(wa.omega3,1)) +wa.alignrV
    h_libero =  -wa.dep*(1.0-min(wa.omega2+wa.omega3,1)) + wa.maxx*(min(wa.omega2+wa.omega3,1)) +wa.alignrV -wa.contiguity -wa.unif*wa.mu
    h_liberalz =  -wa.dep*(1.0-min(wa.beta2,1)) + wa.maxx*(min(wa.beta2,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) +wa.alignrC -wa.complex -wa.nocoda
    h_liberoz = -wa.dep*(1.0-min(wa.omega2+wa.omega3,1)) + wa.maxx*(min(wa.omega2+wa.omega3,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) +wa.alignrC -wa.nocoda -wa.unif*wa.mu -wa.contiguity
    harmonies = [h_liberal, h_liberalo, h_libero, h_liberalz, h_liberoz]
    print h_liberal, h_liberalo, h_libero, h_liberalz, h_liberoz
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 2 in highest:
        print 'Correct result for liberaux', highest
    else:   
        all_correct = False
        print 'Wrong result for liberaux.', highest 
        wa.incomega3(eta1)
        if wa.contiguity > eta1:
            wa.inccontig(-eta1)

    #liberales
    h_liberal =  -wa.dep*(1.0-min(wa.beta2+wa.phi,1)) + wa.maxx*(min(wa.beta2+wa.phi,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) +wa.alignrC -wa.nocoda
    h_liberalo =  -wa.dep*(1.0-min(wa.beta2+wa.phi,1)) + wa.maxx*(min(wa.beta2+wa.phi,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) -wa.dep*(1.0-min(wa.omega3,1)) + wa.maxx*(min(wa.omega3,1)) +wa.alignrV
    h_libero =  -wa.dep*(1.0-min(wa.omega2+wa.omega3+wa.phi,1)) + wa.maxx*(min(wa.omega2+wa.omega3+wa.phi,1)) +wa.alignrV -wa.contiguity -wa.unif*wa.mu
    h_liberalz =  -wa.dep*(1.0-min(wa.beta2+wa.phi,1)) + wa.maxx*(min(wa.beta2+wa.phi,1)) -wa.dep*(1.0-min(wa.omega2,1)) + wa.maxx*(min(wa.omega2,1)) -wa.dep*(1.0-min(wa.zs,1)) + wa.maxx*(min(wa.zs,1)) +wa.alignrC -wa.complex -wa.nocoda
    h_liberoz = -wa.dep*(1.0-min(wa.omega2+wa.omega3,1)) + wa.maxx*(min(wa.omega2+wa.omega3,1)) -wa.dep*(1.0-min(wa.zs+wa.phi,1)) + wa.maxx*(min(wa.zs+wa.phi,1)) +wa.alignrC -wa.nocoda -wa.unif*(1+wa.mu) -wa.contiguity
    harmonies = [h_liberal, h_liberalo, h_libero, h_liberalz, h_liberoz]
    print h_liberal, h_liberalo, h_libero, h_liberalz, h_liberoz
    m = max(harmonies)
    highest = [i for i,j in enumerate(harmonies) if abs(j-m) < 0.001]
    if 0 in highest:
        print 'Correct result for liberales', highest
    else:   
        all_correct = False
        print 'Wrong result for liberales.', highest 
        if wa.omega3 > eta1:
            wa.incomega3(eta1)
        wa.inccontig(eta1)
        wa.incbeta2(eta1)


    print 'beta', wa.beta
    print 'beta2', wa.beta2
    print 'omega', wa.omega
    print 'omega2', wa.omega2
    print 'omega3', wa.omega3
    print 'phi', wa.phi
    print 'contig', wa.contiguity
    print 'complex', wa.complex
    print 'end of loop'
    if all_correct == True:
        print 'all correct'
        break
    else:
        print 'entering new loop'

print 'beta', wa.beta, 'beta2', wa.beta2, 'omega', wa.omega, 'omega2', wa.omega2, 'phi', wa.phi, 'contig', wa.contiguity, 'complex', wa.complex, 'dep', wa.dep, 'maxx', wa.maxx, 't2', wa.t2, 'z2', wa.z2, 'n2', wa.n2, 'combined[0]', wa.combined[0], 'combined[1]', wa.combined[1], 'combined[2]', wa.combined[2], 't1', wa.t1, 'z1', wa.z1, 'zs', wa.zs, 'omega3', wa.omega3, 'n1', wa.n1, 't1f', wa.t1f, 'z1f', wa.z1f, 'n1f', wa.n1f, 'ele', wa.ele, 'eha', wa.eha, 'eho', wa.eho, 'eham', wa.eham, 'ehe', wa.ehe, 'null', wa.null, 'onset', wa.onset, 'unif', wa.unif, 'mumed', mumed, 'mu', wa.mu, 'alignlC', wa.alignlC, 'alignrC', wa.alignrC, 'alignlV', wa.alignlV, 'alignrV', wa.alignrV, 'nocoda', wa.nocoda 

valuedict = {'beta': wa.beta, 'beta2': wa.beta2, 'omega': wa.omega, 'omega2': wa.omega2, 'phi': wa.phi, 'contiguity': wa.contiguity, 'complex': wa.complex, 'dep': wa.dep, 'maxx': wa.maxx, 't2': wa.t2, 'z2': wa.z2, 'n2': wa.n2, 'combined[0]': wa.combined[0], 'combined[1]': wa.combined[1], 'combined[2]': wa.combined[2], 't1': wa.t1, 'z1': wa.z1, 'zs': wa.zs, 'omega3': wa.omega3, 'n1': wa.n1, 't1f': wa.t1f, 'z1f': wa.z1f, 'n1f': wa.n1f, 'ele': wa.ele, 'eha': wa.eha, 'eho': wa.eho, 'eham': wa.eham, 'ehe': wa.ehe, 'null': wa.null, 'onset': wa.onset, 'unif': wa.unif, 'mumed': mumed, 'mu': wa.mu, 'alignlC': wa.alignlC, 'alignrC': wa.alignrC, 'alignlV': wa.alignlV, 'alignrV': wa.alignrV, 'nocoda': wa.nocoda, 'integrity': wa.integrity} 


with open('values.json', 'w') as fp:
    json.dump(valuedict,fp)


print 'finished'
