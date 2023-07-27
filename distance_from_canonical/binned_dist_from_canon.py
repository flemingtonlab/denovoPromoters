import pandas as pd
pwd
ls -tr
upog = pd.read_table('upoG_over3tpmcanonical.bed', header=None,index_col=0)
upog
upog = pd.read_table('upoG_over3tpmcanonical.bed', header=None)
anti = pd.read_table('antisense_over3tpmcanonical.bed', header=None)
anti
%pylab
y = upog[upog[13]<1500]
y
y = upog[upog[13]>-1500]
y
plt.hist(y[13],bins=50,histtype='step')
plt.hist(y[13],bins=500)
upog
upog.sort_values(13)
upog.sort_values(13).iloc[-1]
y
y = upog[(upog[13]>-1500)&(upog[13]<-25)]
y
plt.hist(y[13],bins=500)
plt.hist(y[13],bins=50)
y = upog[(upog[13]>-10000)&(upog[13]<-25)]
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=100)
plt.hist(y[13],bins=200)
plt.hist(y[13],bins=20)
y = upog[(upog[13]>-15000)&(upog[13]<-25)]
plt.hist(y[13],bins=200)
y = upog[(upog[13]>-15000)&(upog[13]<-15)]
plt.hist(y[13],bins=200)
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=70)
plt.hist(y[13],bins=40)
y
y.iloc[3]
y[y[4]>10]
y2=y[y[4]>10]
plt.hist(y2[13],bins=40)
plt.hist(y2[13],bins=40)
plt.hist(y2[13],bins=100)
plt.hist(y2[13],bins=100)
plt.hist(y2[13],bins=40)
plt.hist(y2[13],bins=50)
y2=y[y[4]>20]
plt.hist(y2[13],bins=50)
y2=y[y[4]>8]
plt.hist(y2[13],bins=50)
y2=y[y[4]>=8]
plt.hist(y2[13],bins=50)
y = upog[(upog[13]>-10000)&(upog[13]<-15)]
plt.hist(y2[13],bins=50)
y[4].min
y[4].min()
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=20)
anti
z = anti[(anti[13]<10000)&(anti[13]>15)]
plt.hist(z[13],bins=20)
z
y
z
y
z
y
z
y
y[11].unique()
y[11].value_counts()
plt.hist(z[13],bins=50)
plt.hist(z[13],bins=50)
plt.hist(y[13],bins=50)
pwd
plt.savefig('upog_anti_test.jpg')
pwd
!open .
upog
y
z
y
z
z
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
y
plt.hist(y[13],bins=50)
upog = pd.read_table('upoG_over0tpmcanonical.bed', header=None)
anti = pd.read_table('antisense_over0tpmcanonical.bed', header=None)
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
y
z = anti[(anti[13]<10000)&(anti[13]>1)]
z
y
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(z[13],bins=50)
y
y[12]
y[12].sort_values()
y
int(y[12])
ls -trlh
less upoG_over0tpmcanonical.bed
!less upoG_over0tpmcanonical.bed
ls -tr
x = pd.read_table('TSS_from_mutu_ctl2.bed',sep='\t',header=False)
x = pd.read_table('TSS_from_mutu_ctl2.bed',sep='\t',header=None)
x
x[1]=x[1].map(lambda y:int(y))
x
x[2]=x[2].map(lambda y:int(y))
x[x[6]>5]
x[x[6]>5].to_csv('mutu_canonical.over5tpm.bed',sep='\t')
x[x[6]>5].to_csv('mutu_canonical.over5tpm.bed',sep='\t',header=None)
x[x[6]>5].to_csv('mutu_canonical.over5tpm.bed',sep='\t',header=None,index=None)
ls -trlh
upog = pd.read_table('upoG_over5tpmcanonical.bed', header=None)
anti = pd.read_table('antisense_over5tpmcanonical.bed', header=None)
ls -tr
anti = pd.read_table('antisense_over5tpmcanonical.bed', header=None)
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
z = anti[(anti[13]<10000)&(anti[13]>1)]
plt.hist(z[13],bins=50)
plt.hist(z[13],bins=50)
plt.hist(y[13],bins=50)
y = upog[(upog[13]>-2500)&(upog[13]<-1)]
z = anti[(anti[13]<2500)&(anti[13]>1)]
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(z[13],bins=50)
yu
y
y[12]
len(y[10].unique)
len(y[10].unique())
len(z[10].unique())
y
np.mean(y[12])
y
y[12].map(lambda g:float(g))
y[12]=y[12].map(lambda g:float(g))
z[12]=z[12].map(lambda g:float(g))
np.mean(y[12])
np.mean(z[12])
z
upog
anti
y
y[y[4]>20]
y[y[4]>10]
y
y[y[4]>20]
plt.hist(y[y[4]>20][13],bins=50)
plt.hist(y[y[4]>20][13],bins=50)
plt.hist(z[z[4]>20][13],bins=50)
z
upog
anti
upog
znti
anti
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=20)
plt.hist(y[13],bins=50)
y
upog
upog.groupby(10)
upog.groupby(10).min()[13]
upog.groupby(10).min(13)
upog.groupby(10).max(13)
upog_min = upog.groupby(10).max(13)
y = upog[(upog_min[13]>-2500)&(upog_min[13]<-1)]
y = upog_min[(upog_min[13]>-2500)&(upog_min[13]<-1)]
y
anti_min = upog.groupby(10).min(13)
z = anti_min[(anti_min[13]>-2500)&(anti_min[13]<-1)]
plt.hist(y[13],bins=50)
plt.hist(z[13],bins=50)
z = anti_min[(anti_min[13]<2500)&(anti_min[13]>-1)]
plt.hist(z[13],bins=50)
z
anti_min
anti_min = anti.groupby(10).max(13)
anti_min
plt.hist(z[13],bins=50)
z
z = anti_min[(anti_min[13]<2500)&(anti_min[13]>-1)]
z
plt.hist(z[13],bins=50)
plt.hist(y[13],bins=50)
z = anti_min[(anti_min[13]<10000)&(anti_min[13]>-1)]
y = upog_min[(upog_min[13]>-1000)&(upog_min[13]<-1)]
y = upog_min[(upog_min[13]>-10000)&(upog_min[13]<-1)]
plt.hist(y[13],bins=50)
plt.hist(z[13],bins=50)
y = upog_min[(upog_min[13]>-10000)&(upog_min[13]<-25)]
z = anti_min[(anti_min[13]<10000)&(anti_min[13]>-25)]
plt.hist(z[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(y[13],bins=100)
y
y[13].sort_values()
z[13].sort_values()
y
z
y
plt.hist(y[13],bins=100)
y = upog_min[(upog_min[13]>-10000)&(upog_min[13]<-45)]
plt.hist(y[13],bins=100)
y = upog_min[(upog_min[13]>-10000)&(upog_min[13]<-1)]
z = anti_min[(anti_min[13]<10000)&(anti_min[13]>-1)]
plt.hist(y[13],bins=100)
plt.hist(z[13],bins=100)
plt.hist(z[13],bins=200)
plt.hist(y[13],bins=200)
plt.hist(y[13],bins=200)
plt.hist(y[13],bins=20)
plt.hist(y[13],bins=40)
x
ls -trlh
x
x[x[6]>10]
pwd
x[x[6]>10].to_csv('mutu_canonical_over10tpm.bed',sep='\t',index=None,header=None)
x[x[6]>10].to_csv('mutu_canonical_over10tpm.bed',sep='\t',index=None,header=None)
upog = pd.read_table('upoG_over10tpmcanonical.bed', header=None)
anti = pd.read_table('antisense_over10tpmcanonical.bed', header=None)
anti
upog
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
z = anti[(anti[13]<10000)&(anti[13]>1)]
y
z
plt.hist(y[13],bins=40)
plt.hist(z[13],bins=40)
z
y
y
np.median(y[4])
np.median(z[4])
np.mean(z[4])
np.mean(y[4])
y
y = upog[(upog[13]>-500)&(upog[13]<-1)]
np.mean(y[4])
z
y
z
z = anti[(anti[13]<500)&(anti[13]>1)]
z
y
np.mean(y[4])
np.mean(z[4])
y
y
upog
upog.sort_values(13)
upog.sort_values(13).iloc[-500]
upog.sort_values(13).iloc[-501]
upog.sort_values(13).iloc[-500]
upog = pd.read_table('upoG_over10tpmcanonical.bed', header=None)
upog
anti = pd.read_table('antisense_over10tpmcanonical.bed', header=None)
anti
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
z = anti[(anti[13]<10000)&(anti[13]>1)]
yz
y
z
plt.hist(z[13],bins=40)
plt.hist(z[13],bins=40)
plt.hist(y[13],bins=40)
plt.hist(z[13],bins=200)
plt.hist(y[13],bins=200)
plt.hist(y[13],bins=200)
plt.hist(z[13],bins=200)
plt.hist(y[13],bins=100)
plt.hist(z[13],bins=100)
plt.hist(-z[13],bins=100)
plt.hist(y[13],bins=100)
plt.hist(-z[13],bins=100)
plt.hist(z[13],bins=100)
plt.hist(y[13],bins=100)
x
x[x[6]>=1].to_csv('mutu_canonical_over1tpm.bed',sep='\t',index=None,header=None)
x[x[6]>=1].to_csv('mutu_canonical_over1tpm.bed',sep='\t',index=None,header=None)
upog = pd.read_table('upoG_over1tpmcanonical.bed', header=None)
anti = pd.read_table('antisense_over1tpmcanonical.bed', header=None)
ls -tr
anti = pd.read_table('anti_over1tpmcanonical.bed', header=None)
anti
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
z = anti[(anti[13]<10000)&(anti[13]>1)]
plt.hist(y[13],bins=100)
plt.hist(y[13],bins=100)
y
z
plt.hist(z[13],bins=100)
plt.hist(z[13],bins=50)
plt.hist(y[13],bins=50)
plt.hist(z[13],bins=50)
plt.hist(z[13],bins=100)
plt.hist(y[13],bins=100)
x[x[6]>=5].to_csv('mutu_canonical_over5tpm.bed',sep='\t',index=None,header=None)
upog = pd.read_table('upoG_over5tpmcanonical.bed', header=None)
anti = pd.read_table('anti_over5tpmcanonical.bed', header=None)
y = upog[(upog[13]>-10000)&(upog[13]<-1)]
z = anti[(anti[13]<10000)&(anti[13]>1)]
y
z
y
plt.hist(y[13],bins=100)
plt.hist(y[13],bins=100)
y
z
y
y[y[13]>=-1500]
z[z[13]<=1500]
x
x
x[x[6]>5]
x
y
z
plt.hist(y[13],bins=100)
fig = plt.figure()
ax = plt.subplot()
plt.hist(y[13],bins=100,histtype='step')
plt.hist(z[13],bins=100,histtype='step')
plt.savefig('upog_antisense_near_canonical.svg')
plt.savefig('/Users/nate/upog_antisense_near_canonical.svg')
plt.xlim([-10000,10000])
plt.axvline(-1500,ls='--')
plt.axvline(1500,ls='--')
plt.savefig('/Users/nate/upog_antisense_near_canonical.svg')
plt.hist(z[13],bins=200,histtype='step')
plt.hist(z[13],bins=05,histtype='step')
plt.hist(z[13],bins=50,histtype='step')
plt.close()
z[13].plot.density(color='k', alpha=0.5)
y[13].plot.density(color='k', alpha=0.5)
plt.hist(y[13],bins=100,histtype='step')
plt.xlim([-10000,10000])
plt.hist(y[13],bins=100,histtype='step')
plt.xlim([-10000,10000])
plt.savefig('/Users/nate/upog_near_canonical.svg')
plt.hist(z[13],bins=100,histtype='step')
plt.xlim([-10000,10000])
plt.savefig('/Users/nate/antisense_near_canonical.svg')
plt.hist(y[13],bins=100,histtype='step')
plt.xlim([-10000,10000])
plt.ylim([0,220])
plt.yticks([0,110,220])
plt.savefig('/Users/nate/upog_near_canonical.svg')
plt.hist(z[13],bins=100,histtype='step')
plt.xlim([-10000,10000])
plt.ylim([0,220])
plt.yticks([0,110,220])
plt.savefig('/Users/nate/anti_near_canonical.svg')
upog
history