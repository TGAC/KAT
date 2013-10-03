#!/usr/bin/env python
import sys
from math import *
from scipy.stats import gamma,norm
from scipy import mean,optimize
import numpy as np
def plot_hist(h,points,cap,label=""):
        plot([min(cap,x) for x in h[:points]],label=label)

class KmerSpectra(object):
    """A kmer spectra, comprised of different peaks.
    Contains the general fitting method"""
    ###----------------- INIT, LOAD, DUMP, ETC --------------------
    def __init__(self,histo_file=None,points=10000,column=1,cumulative=False):
        """Init with the objective histogram"""
        self.histogram=[]
        self.peaks=[]
        if histo_file:
            self.read_hist(histo_file,points,column,cumulative=cumulative)

    def read_hist(self,name,points=10000,column=1,cumulative=False):
        f=open(name)
        if cumulative:
            self.histogram=[sum([int(y) for y in x.split()[column:]]) for x in f.readlines() if x and x[0]!='#'][:points][1:]
        else:
            self.histogram=[int(x.split()[column]) for x in f.readlines() if x and x[0]!='#'][:points][1:]
        f.close()
        
    
    def total_values(self,start=1,end=10000):
        return map(sum,zip(*[x.points(start,end) for x in self.peaks]))
    
    ###----------------- FIND AND CREATE PEAKS --------------------
    
    def deriv(self,histogram):
        return [histogram[i+2]-histogram[i] for i in xrange(len(histogram)-2)]
    def progsmoothderiv(self,histogram):
        #always return mean
        return [mean(histogram[i+1:i+1+i/10])-mean(histogram[i-i/10:i+1]) for i in xrange(len(histogram)-2)]
        
    def find_maxima(self, center, radious,min_perc,min_elem,histogram=None):
        if histogram==None:
            hs=self.histogram[center-radious:center+radious]
            #Smoothed df/dx (dx=2 units)
            deriv=self.progsmoothderiv(self.histogram)[center-radious:center+radious]
        else:
            hs=histogram[center-radious:center+radious]
            #Smoothed df/dx (dx=2 units)
            deriv=self.progsmoothderiv(histogram)[center-radious:center+radious]
            
        fmax=hs.index(max(hs))
        if fmax==0 or fmax==len(hs)-1: return 0
        
        #find a single inflection point-> that is the maxima
        #Reject by voting to avoid oscillations.
        failpoints=0
        for i in xrange(fmax):
            if deriv[i]<0: failpoints+=1
        for i in xrange(fmax+1,len(deriv)):
            if deriv[i]>0: failpoints+=1
        if float(failpoints)/(2*radious+1)>.1:
            return 0
        
        #TODO: discard maxima if not 1% of the already contained elements on peaks are there
        
        if sum(hs) < min( (float(min_perc)/100)*sum([x.elements for x in self.peaks]),min_elem):
            print "Distribution on %d too small to be considered (%d elements)" % (center-radious+fmax,sum(hs))
            return 0
        #TODO: do further validation
        #print "maxima found on %d, for %s" % (fmax,hs)
        
        return center-radious+fmax
    
    
    def add_peak_and_update_cuts(self,lm,reset_opt=False):
        #Receives a local maxima, adds a peak there if there was none, updates the cuts.
        #If reset_opt is set, it will reset all the peak counts/sd and only keep the means
        fdists=[x.mean for x in self.peaks]
        for f in fdists:
            if lm>=f-f/5 and lm<=f+f/5:
                print "WARNING!!! Dist on %d is not to be added, because existing dist on %d is too close to it." % (lm,f)
                return False
        fdists.append(lm)
        fdists.sort()
        self.cuts=[self.fmin]
        for i in xrange(len(fdists)-1):
            self.cuts.append(int(fdists[i]*(1+float(fdists[i])/(fdists[i]+fdists[i+1]))))
        self.cuts.append(int(min(len(self.histogram)-1,fdists[-1]*1.5)))
        #print "cuts on %s" % str(self.cuts)
        if reset_opt:
            self.peaks=[]
        for i in xrange(len(fdists)):
            if reset_opt or i==fdists.index(lm):
                self.peaks.append(KmerPeak(lm,fdists[i], fdists[i]/6, sum([self.histogram[j] for j in xrange(self.cuts[i],self.cuts[i+1])])))
                self.peaks.sort(key=lambda x: x.mean)
        return True
    
    def create_peaks(self,min_perc,min_elem):
        fmin=0
        #walk till first local minimum (d(f)>0)
        for i in xrange(len(self.histogram)):
            if self.histogram[i]<self.histogram[i+1]:
                fmin=i
                break
        if not fmin:
            print self.histogram[:25]
            print "Could not find local minima"
            self.fmax=0
            self.maxval=0
            
        #Find the global fm|p(fm)=max(p)
        fmax=self.histogram.index(max(self.histogram[fmin:]))
        #print "fmin=%d (%d) |fmax=%d (%d)" % (fmin,self.histogram[fmin],fmax,self.histogram[fmax])
        if fmax<10:
            print "fmax<10, no analysis can be performed"
            return
        self.fmax=fmax
        self.fmin=fmin
        self.maxval=self.histogram[fmax]
        #Explore fm/x and fm*x for x in [1,4]
        for f in [fmax/4,fmax/3,fmax/2,fmax,fmax*2,fmax*3,fmax*4]:
            if f+f/5<len(self.histogram):
                lm=self.find_maxima(f,f/5,min_perc,min_elem)
                if lm:
                    self.add_peak_and_update_cuts(lm, reset_opt=True)
        #if not fdists:
        #    print "Local maxima not found, skipping further analysis"
        #    return
        #Guess counts for all relevant distributions
        #print "Local maxima found on %s" %(str(fdists))
        
        
    
    ###----------------- OPTIMIZATION --------------------
    
    def optimize_peaks(self,min_perc,min_elem,verbose=False):
        print "------ New independent optimization started"
        sortedpeaks=[x for x in self.peaks]
        sortedpeaks.sort(key=lambda x: -x.elements)
        
        #create a base as the histogram and start from there
        base=[x for x in self.histogram]
        for p in sortedpeaks:
            i=self.peaks.index(p)
            #locally optimize the preak
            print "Optimizing peak: %s on [%d,%d]" % (p,self.cuts[i],self.cuts[i+1])
            if verbose:
                figure()
                plot_hist(base,len(self.histogram),int(self.maxval*1.1),label="Base")
                plot_hist(p.points(1,len(self.histogram)),len(self.histogram),int(self.maxval*1.1),label="before")
            p.constrained_opt(self.cuts[i],base[self.cuts[i]:self.cuts[i+1]])
            if verbose:
                plot_hist(p.points(1,len(self.histogram)),len(self.histogram),int(self.maxval*1.1),label="after")
                legend()
                show()
            print "   Result: %s" % p
            #substract the peak effect from the baseline
            for i in xrange(len(self.histogram)):
                base[i]=base[i]-p.point(i)
            if verbose:
                figure()
                plot_hist(base,len(self.histogram),int(self.maxval*1.1),label="NEW Base")
                legend()
                show()
        updated=False
        for f in [self.fmax/4,self.fmax/3,self.fmax/2,self.fmax,self.fmax*2,self.fmax*3,self.fmax*4]:
            if f+f/5<len(self.histogram) and sum([base[x] for x in xrange(f-f/5,f+f/5)])>.01*sum([p.elements for p in self.peaks]):
                lm=self.find_maxima(f,f/5,min_perc,min_elem,base)
                if lm:
                    if self.add_peak_and_update_cuts(lm,reset_opt=False):
                        updated=True
                #else:
                #    print "No new maxima found on %d +/- %d" % (f,f/5)
        if updated: self.optimize_peaks(min_perc,min_elem,verbose)
        
    def residues(self,p):
        """p has a triplet for each distribution"""
        #print p
        for i in xrange(len(self.peaks)):
            self.peaks[i].mean=p[i*3]
            self.peaks[i].shape=p[i*3+1]
            self.peaks[i].elements=p[i*3+2]
        tv=self.total_values()
        r=[tv[i]-self.histogram[i] for i in xrange(self.cuts[0],self.cuts[-1])]
        #a lot more penalty for being UP:
        for i in xrange(len(r)):
            if r[i]>0:
                r[i]=2*r[i]
        #for j in xrange(len(self.peaks)):
        #    if p[0]>self.peaks[j].target_max*1.1 or p[0]<self.peaks[j].target_max*.9:
        #        #for i in xrange(len(r)): r[i]+=r[i]*10*(float(p[j*3]-self.peaks[j].target_max)/self.peaks[j].target_max)
        #        for i in xrange(len(r)): r[i]+=r[i]*2
        return r
    
    def optimize_overall(self):
        p=[]
        for pk in self.peaks:
            p.append(pk.mean)
            p.append(pk.shape)
            p.append(pk.elements)
        optimize.leastsq(self.residues,p)
        #once the better fit is found, check if by taking the unfitted elements new distributions arise.
        return
                   
                   
        
class KmerPeak(object):
    """A distribution representing kmers covered a certain number of times.
    Contains methods for fitting to an interval"""
    
    def __init__(self,target_max,mean,shape,elements):
        self.target_max=target_max
        self.mean=mean
        self.shape=shape
        
        
        self.elements=elements
        
    def __str__(self):
        return "<Kmer distribution peak on %d, with %d elements>" % (self.mean,self.elements)
        
    def point(self,x):
        """Normalized Gaussian"""
        return float(self.elements) / np.sqrt(2 * np.pi) / self.shape * np.exp(-(x - self.mean) ** 2 / 2. / self.shape ** 2)
    
    def points(self,start,end):
        return [self.point(x) for x in xrange(start,end)]
    
    def residues(self,p,offset,histogram):
        """residues of fitting self with parameters p to offset obj[0] with values obj[1]"""
        #TODO: enforce a constrain?
        self.mean=p[0]
        self.shape=p[1]
        self.elements=p[2]
        r=[self.point(offset+i)-histogram[i] for i in xrange(len(histogram))]
        #a lot more penalty for being UP:
        for i in xrange(len(r)):
            if r[i]>0:
                r[i]=2*r[i]
        #if p[0]>self.target_max*1.1 or p[0]<self.target_max*.9:
            #for i in xrange(len(r)): r[i]+=r[i]*10*(float(p[0]-self.target_max)/self.target_max)
        #    for i in xrange(len(r)): r[i]+=r[i]*2
        #penalizatin for moving the mean
        
        return r
    
    def constrained_opt(self,offset,histogram):
        """Does a constrained optimization of fitting to match the histogram"""
        #print "Optimizing peak %s" % str(self)
        optimize.leastsq(self.residues,(self.mean,self.shape,self.elements),(offset,histogram))
        #self.mean=p[0]
        #self.shape=p[1]
        #self.elements=p[2]
        pass
    

class KmerSpectraAnalysis(object):
    
    def __init__(self,filename,points=10000):
        self.spectra=KmerSpectra(filename,points)
        
    def plot_hist(self,h,points,cap):
        plot([min(cap,x) for x in h[:points]])
    
    def plot_hist_df(self,h,points,cap):
        plot([max(-cap,min(cap,x)) for x in [h[i+1]-h[i] for i in xrange(points)]])
    
    def plot_all(self,points,cap):
        self.plot_hist(self.spectra.histogram,points,cap)
        
        #self.plot_hist_df(self.histogram,points,cap)
        self.plot_hist(self.spectra.total_values(1,points+1),points,cap)
    
    def plot_peaks(self,points,cap):
        for p in self.spectra.peaks:
            self.plot_hist(p.points(1,points+1),points,cap)
    
    def analyse(self):
        self.spectra.create_peaks()
        limy=int(self.spectra.maxval*1.1/1000)*1000
        limx=self.spectra.peaks[-1].mean*2
        print "Plot limits: y->%d, x->%d" % (limy,limx)
        self.plot_all(limx,limy)
        show()
        self.spectra.optimize_peaks()
        figure()
        self.plot_all(limx,limy)
        show()
        self.spectra.optimize_overall()
        figure()
        self.plot_all(limx,limy)
        figure()
        self.plot_peaks(limx,limy)
        show()
        for p in self.spectra.peaks: print p

class MXKmerSpectraAnalysis(object):
    def __init__(self,filename,columns=3,points=10000):
        self.spectras=[KmerSpectra(filename,points,column=0,cumulative=True)]
        for i in xrange(columns):
            self.spectras.append(KmerSpectra(filename,points,column=i,cumulative=(i==columns-1)))
        
        
    def plot_hist(self,h,points,cap,label=""):
        plot([min(cap,x) for x in h[:points]],label=label)
        xlabel('Kmer Frecuency')
        ylabel('Kmer Count')
        legend()
    
    def plot_hist_df(self,h,points,cap):
        plot([max(-cap,min(cap,x)) for x in [h[i+1]-h[i] for i in xrange(points)]])
    
    def plot_all(self,points=0,cap=0,spectra=True,fit=True,dists=True):
        if 0==points: points=self.limx
        if 0==cap: cap=self.limy
        for s in self.spectras:
            if self.spectras[0]==s:
                slabel="General Spectra"
            else:
                slabel="%d x present" % (self.spectras.index(s)-1)
            figure()
            if spectra: self.plot_hist(s.histogram,points,cap,label=slabel)
            if fit: self.plot_hist(s.total_values(1,points+1),points,cap,label=slabel+" fit")
            if dists:
                for p in s.peaks:
                    self.plot_hist(p.points(1,points+1),points,cap,label="fit dist %d" % s.peaks.index(p))
            show()
            for p in s.peaks:
                print p
    
    def analyse(self,min_perc=1,min_elem=100000,verbose=False):
        self.limx=0
        self.limy=0
        for s in self.spectras:
            print "analysing spectra... ",
            sys.stdout.flush()
            s.create_peaks(min_perc=min_perc,min_elem=min_elem)
            
            sys.stdout.flush()
            if s.peaks:
                self.limy=max(int(s.maxval*1.1/1000)*1000, self.limy )
                self.limx=max(min(s.peaks[-1].mean*2,len(s.histogram)) , self.limx)
                print "peaks created ... ",
                sys.stdout.flush()
                s.optimize_peaks(min_perc=min_perc,min_elem=min_elem,verbose=verbose)
                print "locally optimised ... ",
                for p in s.peaks: print p
                sys.stdout.flush()
                s.optimize_overall()
                print "overall optimised ... DONE"
                sys.stdout.flush()
        print "Plot limits: y->%d, x->%d" % (self.limy,self.limx)
        
    def peak_stats(self):
        """TODO: Runs analyse (TODO:include general spectra)
                 Takes enough peaks as to cover a given % of the elements:
                     - Find the peak across all distributions
                     - Reports peak stats
                 If multiple peaks have been analyzed, tries to find the "main unique" and explains the results based on that freq.
        """
        #step 1, try to find a reasonable mean for kmer frequency.
        #weighted means by number of elements?
        general_dists=[(x.mean,x.elements) for x in self.spectras[0].peaks]
        general_dists.sort(key=lambda x: x[0])
        pkc=general_dists[0][0]
        kcov=sum([ (x[0]/round(x[0]/pkc) )*x[1] for x in general_dists])/sum([x[1] for x in general_dists])
        print pkc, kcov
        #step 2, selects frequencies for peaks from bigger to smaller till X% of the elements are covered or no more peaks
        goal=0.99*sum([x[1] for x in general_dists])
        maxpeaks=10
        general_dists.sort(key=lambda x: -x[1])
        af=[]
        peaks=0
        covered=0
        for x in general_dists:
            af.append(kcov*round(x[0]/kcov))
            peaks+=1
            covered+=x[1]
            if peaks==maxpeaks or covered>goal:
                break
        #step 3, report for each peak
        #get the candidate peak on each spectra
        for f in af:
            total=0
            pd={}
            for i in xrange(len(self.spectras)-1):
                m=[(x.mean,x.elements) for x in self.spectras[1+i].peaks if x.mean>0.8*f and x.mean<1.2*f]
                if len(m)==1:
                    pd[i]=m[0]
                    total+=m[0][1]
                if len(m)>1:
                    print "WARNING, MORE THAT 1 PEAK FOR f=%.3f FOUND ON THE %dx SPECTRA!!!" % (f,i)
            print "\n---- Report for f=%.3f (total elements %d)----" % (f,total)
            for i in xrange(len(self.spectras)-1):
                if i in pd.keys():
                    print " %dx: %.2f%% (%d elements on f=%.2f)" % (i,float(pd[i][1])*100/total,pd[i][1],pd[i][0])
                else:
                    print " %dx: No significant content" % i
            
        #step 4, general report
        return

if __name__ == '__main__':
    if len(sys.argv)!=6:
        print "Usage: %s filename distribution_count fmax min_perc min_elem" % sys.argv[0]
        exit(1)

    a=MXKmerSpectraAnalysis(sys.argv[1],int(sys.argv[2]),int(sys.argv[3]))
    a.analyse(min_perc=int(sys.argv[4]),min_elem=int(sys.argv[5]),verbose=False)
    a.peak_stats() 
