import numpy as np
from scipy.stats import pearsonr
def returnSNPvals(size1,freq):
    scoresA1=np.random.choice([0,1],size=size1,p=[1-freq,freq])
    scoresA2=np.random.choice([0,1],size=size1,p=[1-freq,freq])
    scoresA=scoresA1+scoresA2
    return scoresA
def returnSNPvalsIfDeviate(size1,fracs):
    num0s,num1s=int(round((size1*fracs[0]))),int(round((size1*fracs[1])))
    num2s=size1-num0s-num1s
    scores0,scores1=np.zeros((num0s),int),np.ones((num1s),int)
    scores2=2*np.ones((num2s),int)
    scoresA1=np.concatenate((scores0,scores1,scores2))
    perm=np.random.permutation(size1)
    scoresA=scoresA1[perm]
    return scoresA
def returnNorm(listScores,varyScores,idx,corrsAll):
    listScoresTemp=listScores.copy()
    listScoresTemp[idx]=varyScores
    normTemp=0
    normStore=[]
    num0=len(corrsAll)
    for key,corr0 in corrsAll.items(): 
        i,j=key
        corrTemp,_= pearsonr(listScoresTemp[i],listScoresTemp[j])
        diff=corr0-corrTemp
        normVal=np.linalg.norm(diff)
        normTemp=normTemp+normVal
        normStore.append(normVal)
    normTemp=normTemp/num0
    return normTemp,normStore
def checkUpDate(scoresI,scoresList,corrsAll,i,norm0):
    num0=int(scoresI.shape[0]/2)
    randNums=np.random.permutation(2*num0)
    for k in range(num0):
        idx1,idx2=randNums[2*k],randNums[2*k+1]
        score1,score2=scoresI[idx2],scoresI[idx1]
        if score1-score2!=0:
            scoresTemp=scoresI.copy()
            scoresTemp[idx1],scoresTemp[idx2]=score1,score2
            normTemp,normStore=returnNorm(scoresList,scoresTemp,i,corrsAll)
            if normTemp<norm0:
                scoresI=scoresTemp
                norm0=normTemp
    return norm0,scoresI

def returnUpdatedScores(scoresList,size1,corrsAll,tolerance):
    keepRunning=True
    # tolerance=0.01
    listLen=len(scoresList)
    # print('Length of list=',listLen)
    iter0=0
    # print('correlations to match=',corrsAll)
    rand0=np.random.permutation(listLen)
    averageNums=[]
    for i in range(listLen):
        averageNums.append(np.mean(scoresList[i]))
    # print(averageNums)
    
    
    while keepRunning:
        for i0 in range(listLen): 
            # print(listLen)
            i=rand0[i0]
            scoresI=scoresList[i] 
            norm0,normStore=returnNorm(scoresList,scoresI,i,corrsAll)
            norm0,scoresI=checkUpDate(scoresI,scoresList,corrsAll,i,norm0)
            scoresList[i]=scoresI
            if norm0<tolerance:
                keepRunning=False
                success=True
        # print(normStore)
        iter0+=1
        if iter0%200==0:
            randInt,randomise=np.random.randint(listLen),np.random.permutation(size1)
            scoresRand=scoresList[randInt][randomise]
            scoresList[randInt]=scoresRand
            # print(iter0)
        if iter0==10000:
            keepRunning=False
            success=False       
    return scoresList,success














