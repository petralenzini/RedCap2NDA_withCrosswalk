
#two structures, ndarsubjects01 and edinburgh_hand01 are created by specialty
# programs HCD_ndar_edinburgh.py
#each of these structures requires consideration of all HCPD-redcap databases
#and were used as guinea pig behavioral data structures
import os, datetime
import json
import sys
from multiprocessing.dummy import Pool
import pandas as pd
import pycurl
import numpy as np
import io
from io import BytesIO
import download.box
from download.box import LifespanBox
box_temp='/home/petra/UbWinSharedSpace1/boxtemp'
box = LifespanBox(cache=box_temp)
#cache_space='/home/petra/UbWinSharedSpace1/scratch'
cache_space=box_temp

from boxsdk import JWTAuth, OAuth2, Client
snapshotdate = datetime.datetime.today().strftime('%m_%d_%Y')

#cur_dir = os.path.dirname(os.path.abspath(__file__))
#default_cache = "/data/intradb/tmp/box2nda_cache"

###REDCAP API tokens moved to configuration file
redcapconfigfile="/home/petra/UbWinSharedSpace1/ccf-nda-behavioral/PycharmToolbox/.boxApp/redcapconfig.csv"
#redcapconfigfile='/data/intradb/home/.boxApp/redcapconfig.csv'
boxconfigfile="/home/petra/UbWinSharedSpace1/ccf-nda-behavioral/PycharmToolbox/.boxApp/config.json"
crosswalkfile="/home/petra/UbWinSharedSpace1/redcap2nda_Lifespan2019/HCA_crosswalk_docs/HCPA_Crosswalk_concatenated_25Nov2019.csv"
#crosswalkfile="/home/petra/UbWinSharedSpace1/redcap2nda_Lifespan2019/crosswalk_docs/Crosswalk_test.csv"
crosswalk=pd.read_csv(crosswalkfile)

ndar_fields='UnrelatedHCAHCD_w_STG_Image_and_pseudo_GUID09_27_2019.csv'
ndar=pd.read_csv('/home/petra/UbWinSharedSpace1/redcap2nda_Lifespan2019/Dev_pedigrees/'+ndar_fields)
ndar=ndar.loc[ndar.subjectped.str.contains('HCA')]
pathout="/home/petra/UbWinSharedSpace1/redcap2nda_Lifespan2019/HCA_crosswalk_docs/prepped_structures"

#will import the crosswalk and sort by structures, then do for all structures the following
#get fieldlist for items in this structure (if redcap based and not special case like session drugs)
#get the data from redcap corresponding to this fieldlist (sincerely hope no merging required with Box data...
#needs to be called study data because code snpts in crosswalk refer to studydata datafram
#structures=crosswalk.groupby('NDA Structure').count().reset_index()[['NDA Structure']]
structures=crosswalk.drop_duplicates(subset='NDA Structure')[['HCP-A Source','dbase','NDA Structure','specialty_code']]

#normal redcap structures
normals=structures.loc[(structures['HCP-A Source'].str.contains('Redcap')==True) & (structures.specialty_code.isnull()==True)]

for structure in normals['NDA Structure']:
    elements=crosswalk.loc[crosswalk['NDA Structure']==structure][['HCP-A Element']]
    vars=list(elements['HCP-A Element'])
    study=crosswalk.loc[crosswalk['NDA Structure']==structure,'dbase'].values[0]
    redcap2structure(vars,crosswalk,pathstructuresout=pathout,studystr=study)



#ravlt - special case
extras=getredcapfieldsjson(fieldlist=[],study='hcpa')
crosswalk_subset=crosswalk.loc[(crosswalk['HCP-A Source'].str.contains('Box'))
                    & (crosswalk['NDA Structure']=='ravlt01')]

ravltid=crosswalk_subset.dbase.unique()[0]
ravlt=Box2dataframe(ravltid)
ravlt=ravlt.loc[ravlt.visit=='V1']
#execute any specialty codes
for index,row in crosswalk_subset.iterrows():
    if pd.isna(row['python first code for form request']):
        pass
    else:
        exec(row['python first code for form request'])
listvars=crosswalk_subset['HCP-A Element name in uploaded file'].tolist()
ravlt=ravlt[['subject']+listvars]
studydata=pd.merge(ravlt,extras,on='subject',how='right')
transformed=studydata.loc[studydata.flagged.isnull()==True].drop(columns='flagged')
#merge with required fields from ndar subjects
#--age and date need be recalcuated on fly from redcap data because possiblity of multiple dates per person'
ndarsub=ndar[['nda_guid','subjectped']].rename(
    columns={'nda_guid':'subjectkey','subjectped':'src_subject_id'}).copy()
dout=pd.merge(ndarsub,transformed,how='left',left_on='src_subject_id',right_on='subject').drop(columns={'subject','dob','site', 'study', 'subject_id'})
crosswalk_subset.reset_index(inplace=True)
strucroot=crosswalk_subset['NDA Structure'].str.strip().str[:-2][0]
strucnum=crosswalk_subset['NDA Structure'].str.strip().str[-2:][0]
#finalsubset - i.e. no withdraws
#subjectkey	src_subject_id	interview_age	interview_date	gender
filePath=os.path.join(pathout,'HCPA_'+strucroot+strucnum+'_'+snapshotdate+'.csv')
if os.path.exists(filePath):
    os.remove(filePath)
else:
    pass
    #print("Can not delete the file as it doesn't exists")

with open(filePath,'a') as f:
    f.write(strucroot+","+str(int(strucnum))+"\n")
    dout.to_csv(f,index=False)


#######################################
### special cases  ###
#########################################
#medhx vars merge from hcpa and ssaga databases
crosswalk_subset=crosswalk.loc[crosswalk['NDA Structure']=='medh01']
hcpalist=crosswalk_subset.loc[crosswalk_subset.dbase=='hcpa','HCP-A Element'].tolist()
hcpadata=getredcapfieldsjson(fieldlist=hcpalist, study='hcpa')
ssagalist=crosswalk_subset.loc[crosswalk_subset.dbase=='ssaga','HCP-A Element'].tolist()
ssagadata=getredcapfieldsjson(fieldlist=ssagalist, study='ssaga')
ssagahcpa=pd.merge(hcpadata,ssagadata.drop(columns=['flagged','study','interview_date']),on='subject',how='inner')
vars=hcpalist+ssagalist
redcap2structure(vars,crosswalk_subset,pathstructuresout=pathout,studystr='hcpa',dframe=ssagahcpa)




#caffeine nicotine and other drug sessions
#need six rows per person corresponding to 6 sessions
crosswalk_subset=crosswalk.loc[crosswalk['NDA Structure']=='drugscr01']
sessions=['1','2','3','4','5','6']
renamelist=crosswalk_subset.loc[(crosswalk_subset['HCP-A Element'].str.contains('s1')==True) |
    (crosswalk_subset['HCP-A Element'].str.contains('drug1')==True),'HCP-A Element'].tolist() + ['alc_breath1']
allsessions=pd.DataFrame()
for session in sessions:
    slist=crosswalk_subset.loc[(crosswalk_subset['HCP-A Element'].str.contains('s'+session)==True) |
        (crosswalk_subset['HCP-A Element'].str.contains('drug'+session)==True),'HCP-A Element'].tolist() + ['alc_breath'+session]
    allstudydata=pd.DataFrame()
    studies = crosswalk_subset.dbase.values[0]
    for pop in studies.split():
        studydata = pd.DataFrame()
        studydata = getredcapfieldsjson(fieldlist=slist, study=pop)
        allstudydata = pd.concat([allstudydata, studydata], axis=0, sort=True)
    allstudydata['caffeine_s' + session + '___1'] = allstudydata['caffeine_s' + session + '___1'].str.replace('1','Prior to visit')
    allstudydata['caffeine_s' + session + '___2'] = allstudydata['caffeine_s' + session + '___2'].str.replace('2','During visit')
    allstudydata['caffeine_s' + session]=allstudydata['caffeine_s'+session+'___1'] + ' ' + allstudydata['caffeine_s'+session+'___2']
    allstudydata['caffeine_s' + session]=allstudydata['caffeine_s' + session].str.replace('0','')
    allstudydata['nicotine_s' + session + '___1'] = allstudydata['nicotine_s' + session + '___1'].str.replace('1','Prior to visit')
    allstudydata['nicotine_s' + session + '___2'] = allstudydata['nicotine_s' + session + '___2'].str.replace('2','During visit')
    allstudydata['nicotine_s' + session]=allstudydata['nicotine_s'+session+'___1'] + ' ' + allstudydata['nicotine_s'+session+'___2']
    allstudydata['nicotine_s' + session]=allstudydata['nicotine_s' + session].str.replace('0','')
    allstudydata=allstudydata.drop(columns=['caffeine_s'+session+'___1','caffeine_s'+session+'___2','nicotine_s'+session+'___1','nicotine_s'+session+'___2']).copy()
    varmap = dict(zip(slist, renamelist))
    allstudydata=allstudydata.rename(columns=varmap)
    allstudydata['version_form']=session
    allsessions=pd.concat([allsessions,allstudydata],axis=0,sort=True)

lout=list(crosswalk_subset['HCP-A Element name in uploaded file'])
cleanedlist = [x for x in lout if str(x) != 'nan']
listout=['subject','flagged','interview_date','interview_age','gender']+cleanedlist  #output new variables and subset to those not flagged for withdrawal.
transformed=allsessions[listout].loc[allsessions.flagged.isnull()==True].drop(columns={'flagged','interview_date','gender','interview_age'})
#merge with required fields from vars in intradb staging (guid, etc)
#not sure whether it makes sense to pull these in here or recalculate on fly from redcap.
#future issues:  compare this approach (e.g. pull from the file above named 'ndar') vs. what happens in the applycrosswalk.py
#program for HCD, which regenerates on fly...will require some recodeing below to pull from redcap...
#might just be easier to pull once...but how will this affect visit numbers?
ndarsub=ndar[['nda_guid','subjectped','nda_gender','nda_interview_age','nda_interview_date']].rename(
    columns={'nda_guid':'subjectkey','subjectped':'src_subject_id','nda_gender':'gender',
             'nda_interview_date':'interview_date','nda_interview_age':'interview_age'}).copy()
dout=pd.merge(ndarsub,transformed,how='left',left_on='src_subject_id',right_on='subject').drop(columns='subject')
dout['interview_date'] = pd.to_datetime(dout['interview_date']).dt.strftime('%m/%d/%Y')
#now export
crosswalk_subset.reset_index(inplace=True)
strucroot=crosswalk_subset['NDA Structure'].str.strip().str[:-2][0]
strucnum=crosswalk_subset['NDA Structure'].str.strip().str[-2:][0]
#finalsubset - i.e. no withdraws
#subjectkey	src_subject_id	interview_age	interview_date	gender
filePath=os.path.join(pathout,'HCPA_'+strucroot+strucnum+'_'+snapshotdate+'.csv')
if os.path.exists(filePath):
    os.remove(filePath)
else:
    pass
    #print("Can not delete the file as it doesn't exists")

    with open(filePath,'a') as f:
        f.write(strucroot+","+str(int(strucnum))+"\n")
        dout.to_csv(f,index=False)
######################################

#create the two penncnp stuctures
structuresbox=crosswalk_subset.drop_duplicates(subset='NDA Structure')[['HCP-A Source','dbase','NDA Structure','specialty_code']]
for structure in structuresbox['NDA Structure']:
    listvars=crosswalk_subset.loc[crosswalk_subset['NDA Structure']==structure,'HCP-A Element name in uploaded file'].tolist()
    pennstruct=penn[['subid']+listvars]
    studydata=pd.merge(pennstruct,extras,left_on='subid',right_on='subject',how='right')
    transformed=studydata.loc[studydata.flagged.isnull()==True].drop(columns='flagged')
    #merge with required fields from ndar subjects
    #--age and date need be recalcuated on fly from redcap data because possiblity of multiple dates per person'
    ndarsub=ndar[['nda_guid','subjectped']].rename(
    columns={'nda_guid':'subjectkey','subjectped':'src_subject_id'}).copy()
    dout=pd.merge(ndarsub,transformed,how='left',left_on='src_subject_id',right_on='subject').drop(columns={'subject','dob','site', 'study', 'subject_id'})
    crosswalk_boxsubset=crosswalk_subset.loc[crosswalk_subset['NDA Structure']==structure]
    crosswalk_boxsubset.reset_index(inplace=True)
    strucroot=crosswalk_boxsubset['NDA Structure'].str.strip().str[:-2][0]
    strucnum=crosswalk_boxsubset['NDA Structure'].str.strip().str[-2:][0]
    filePath=os.path.join(pathout,'HCPA_'+strucroot+strucnum+'_'+snapshotdate+'.csv')
    if os.path.exists(filePath):
        os.remove(filePath)
    else:
        pass
        #print("Can not delete the file as it doesn't exists")
    with open(filePath,'a') as f:
        f.write(strucroot+","+str(int(strucnum))+"\n")
        dout.to_csv(f,index=False)


###end penncnp
######################################


def Box2dataframe(curated_fileid_start):#,study,site,datatype,boxsnapshotfolderid,boxsnapshotQCfolderid):
    #get current best curated data from BOX and read into pandas dataframe for QC
    raw_fileid=curated_fileid_start
    rawobject=box.download_file(raw_fileid)
    data_path = os.path.join(cache_space, rawobject.get().name)
    raw=pd.read_csv(data_path,header=0,low_memory=False, encoding = 'ISO-8859-1')
    #raw['DateCreatedDatetime']=pd.to_datetime(raw.DateCreated).dt.round('min')
    #raw['InstStartedDatetime']=pd.to_datetime(raw.InstStarted).dt.round('min')
    #raw['InstEndedDatetime']=pd.to_datetime(raw.InstEnded).dt.round('min')
    return raw



#use json format because otherwise commas in strings convert wrong in csv read
def getredcapfieldsjson(fieldlist, study='hcpa'):  # , token=token[0],field=field[0],event=event[0]):
    """
    Downloads requested fields from Redcap databases specified by details in redcapconfig file
    Returns panda dataframe with fields 'study', 'Subject_ID, 'subject', and 'flagged', where 'Subject_ID' is the
    patient id in the database of interest (sometimes called subject_id, parent_id) as well as requested fields.
    subject is this same id stripped of underscores or flags like 'excluded' to make it easier to merge
    flagged contains the extra characters other than the id so you can keep track of who should NOT be uploaded to NDA
    or elsewwhere shared
    """
    auth = pd.read_csv(redcapconfigfile)
    studydata = pd.DataFrame()
    fieldlistlabel = ['fields[' + str(i) + ']' for i in range(5, len(fieldlist) + 5)]
    fieldrow = dict(zip(fieldlistlabel, fieldlist))
    d1 = {'token': auth.loc[auth.study == study, 'token'].values[0], 'content': 'record', 'format': 'json', 'type': 'flat',
          'fields[0]': auth.loc[auth.study == study, 'field'].values[0],
          'fields[1]': auth.loc[auth.study == study, 'interview_date'].values[0],
          'fields[2]': auth.loc[auth.study == study, 'sexatbirth'].values[0],
          'fields[3]': auth.loc[auth.study == study, 'sitenum'].values[0],
          'fields[4]': auth.loc[auth.study == study, 'dobvar'].values[0]}
    d2 = fieldrow
    d3 = {'events[0]': auth.loc[auth.study == study, 'event'].values[0], 'rawOrLabel': 'raw', 'rawOrLabelHeaders': 'raw',
          'exportCheckboxLabel': 'false',
          'exportSurveyFields': 'false', 'exportDataAccessGroups': 'false', 'returnFormat': 'json'}
    data = {**d1, **d2, **d3}
    buf = BytesIO()
    ch = pycurl.Curl()
    ch.setopt(ch.URL, 'https://redcap.wustl.edu/redcap/srvrs/prod_v3_1_0_001/redcap/api/')
    ch.setopt(ch.HTTPPOST, list(data.items()))
    ch.setopt(ch.WRITEDATA, buf)
    ch.perform()
    ch.close()
    htmlString = buf.getvalue().decode('UTF-8')
    buf.close()
    d = json.loads(htmlString)
    #parent_ids = pd.DataFrame(htmlString.splitlines(), columns=['row'])
    #header = parent_ids.iloc[0]
    #headerv2 = header.str.replace(auth.loc[auth.study == study, 'interview_date'].values[0], 'interview_date')
    #headerv3 = headerv2.str.split(',')
    #parent_ids.drop([0], inplace=True)
    #pexpanded = pd.DataFrame(parent_ids.row.str.split(pat='\t').values.tolist(), columns=headerv3.values.tolist()[0])
    pexpanded=pd.DataFrame(d)
    pexpanded = pexpanded.loc[~(pexpanded[auth.loc[auth.study == study, 'field'].values[0]] == '')]  ##
    new = pexpanded[auth.loc[auth.study == study, 'field'].values[0]].str.split("_", 1, expand=True)
    pexpanded['subject'] = new[0].str.strip()
    pexpanded['flagged'] = new[1].str.strip()
    pexpanded['study'] = study  # auth.study[i]
    studydata = pd.concat([studydata, pexpanded], axis=0, sort=True)
    studydata=studydata.rename(columns={auth.loc[auth.study == study, 'interview_date'].values[0]:'interview_date'})
    # Convert age in years to age in months
    # note that dob is hardcoded var name here because all redcap databases use same variable name...sue me
    # interview date, which was originally v1_date for hcpa, has been renamed in line above, headerv2
    try:
        studydata['nb_months'] = (
                12 * (pd.to_datetime(studydata['interview_date']).dt.year - pd.to_datetime(studydata.dob).dt.year) +
                (pd.to_datetime(studydata['interview_date']).dt.month - pd.to_datetime(studydata.dob).dt.month) +
                (pd.to_datetime(studydata['interview_date']).dt.day - pd.to_datetime(studydata.dob).dt.day) / 31)
        studydatasub=studydata.loc[studydata.nb_months.isnull()].copy()
        studydatasuper = studydata.loc[~(studydata.nb_months.isnull())].copy()
        studydatasuper['nb_months'] = studydatasuper['nb_months'].apply(np.floor).astype(int)
        studydatasuper['nb_monthsPHI'] = studydatasuper['nb_months']
        studydatasuper.loc[studydatasuper.nb_months > 1080, 'nb_monthsPHI'] = 1200
        studydata=pd.concat([studydatasub,studydatasuper],sort=True)
        studydata = studydata.drop(columns={'nb_months'}).rename(columns={'nb_monthsPHI': 'interview_age'})
    except:
        pass
    #convert gender to M/F string
    try:
        studydata.gender = studydata.gender.str.replace('1', 'M')
        studydata.gender = studydata.gender.str.replace('2', 'F')
    except:
        print(study+' has no variable named gender')
    return studydata


def redcap2structure(vars,crosswalk,pathstructuresout=pathout,studystr='hcpa',dframe=None):
    """
    Takes list of vars from the crosswalk, gets the data from Redcap, and puts into structure format after
    merging with NDAR requiredvars.  Outputs a csv structure in NDA format to pathstructureout location
    """
    #varslim=[x for x in fieldlist if str(x) != 'nan']
    #varnan=[x for x in fieldlist if str(x) == 'nan']
    if dframe is not None:
        studydata=dframe
    else:
        studydata=getredcapfieldsjson(fieldlist=vars,study=studystr)
    #studydata=getredcapfieldsjson(fieldlist=vars,study=studystr)
    if studystr=='ssaga':
        extras=getredcapfieldsjson(fieldlist=[],study='hcpa')
        studydata=pd.merge(studydata.drop(columns='interview_date'),extras[['interview_age','interview_date','subject','gender']],on='subject',how='left')
    #get the relevant rows of the crosswalk
    crosswalk_subset=pd.merge(crosswalk,pd.DataFrame(vars,columns=['HCP-A Element']),on='HCP-A Element',how='inner')[['NDA Structure', 'NDA Element', 'HCP-A Element', 'HCP-A Source',
       'CCF action applied (e.g. request from Form Request)',
       'HCP-A Element name in uploaded file',
       'python first code for form request']]
    #execute transformation codes stored in the crosswalk
    for index,row in crosswalk_subset.iterrows():
        if pd.isna(row['python first code for form request']):
            pass
        else:
            exec(row['python first code for form request'])
    #remove fields with empty values HCP-A Element name in uploaded file -- these are empty because NDA doesnt want them
    crosswalk_subset=crosswalk_subset.loc[crosswalk_subset['HCP-A Element name in uploaded file'].isnull()==False]
    listout=['subject','flagged','interview_date','interview_age','gender']+list(crosswalk_subset['HCP-A Element name in uploaded file'])
    #output new variables and subset to those not flagged for withdrawal.
    transformed=studydata[listout].loc[studydata.flagged.isnull()==True].drop(columns='flagged')
    #merge with required fields from ndar subjects
    #--age and date need be recalcuated on fly from redcap data because possiblity of multiple dates per person'
    ndarsub=ndar[['nda_guid','subjectped']].rename(
        columns={'nda_guid':'subjectkey','subjectped':'src_subject_id'}).copy()
    dout=pd.merge(ndarsub,transformed,how='left',left_on='src_subject_id',right_on='subject').drop(columns='subject')
    dout['interview_date'] = pd.to_datetime(dout['interview_date']).dt.strftime('%m/%d/%Y')
    #now export
    crosswalk_subset.reset_index(inplace=True)
    strucroot=crosswalk_subset['NDA Structure'].str.strip().str[:-2][0]
    strucnum=crosswalk_subset['NDA Structure'].str.strip().str[-2:][0]
    #finalsubset - i.e. no withdraws
    #subjectkey	src_subject_id	interview_age	interview_date	gender
    filePath=os.path.join(pathstructuresout,'HCPA_'+strucroot+strucnum+'_'+snapshotdate+'.csv')
    if os.path.exists(filePath):
        os.remove(filePath)
    else:
        pass
        #print("Can not delete the file as it doesn't exists")

    with open(filePath,'a') as f:
        f.write(strucroot+","+str(int(strucnum))+"\n")
        dout.to_csv(f,index=False)

