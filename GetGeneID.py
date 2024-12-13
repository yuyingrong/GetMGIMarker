# YR 20220509
# Python 3.8+
# this script scrapes Rhesus monkey gene ID and Ensembl ID from NCBI Gene website
# (then converts gene ID to Ensembl ID) -- commented out this step because this could be done in R with package "biomaRt"


from selenium import webdriver
from selenium.webdriver.chrome.service import Service
#from time import sleep


# install chromedriver on MacOS M1
# followed: https://www.tutorialspoint.com/how-to-setup-chrome-driver-with-selenium-on-macos
# ChromeDriver downloaded from: https://chromedriver.chromium.org/downloads

# path to chromedriver installed on my machine
chromedriver_path = '/usr/local/bin/chromedriver'

# referenced: https://blog.csdn.net/m0_62298204/article/details/120802053
# specify path to web driver/browser
driver = webdriver.Chrome(service=Service(chromedriver_path))
# warning if use executable_path=chromedriver_path, which is deprecated
# instead, use service=turn the above path into a Service Object with Service()


# output dict: list(monkey Gene ID, monkey Ensembl ID)
d={}

# input list: human markers
l=['CD25','CD32','CD64','CD86','CD127','CD215','IRF5','iNO5']

def convert(InList, OutDict, organism):

    for name in InList:
        # open a Chrome browser to intended page
        driver.get('https://www.ncbi.nlm.nih.gov/gene')
        #driver.implicitly_wait(10) - no longer needed bc now driver will wait until loaded
        #print(driver.title)
    
        # create a tuple to store output gene ID and Ensembl ID
        OutDict[name]={}

        # locate the textarea, whose <id> is 'term'
        # elem=driver.find_element_by_id('term')
        # the above form is deprecated; use form below
        elem=driver.find_element(by='id',value='term')
        # clear the previous search term if any
        elem.clear()
        
        # type protein name into textarea
        elem.send_keys(f'{name} {organism}')
        # submit the form containing this <id>; i.e. search the protein ID
        elem.submit()

        # locate the first search result
        # elem=driver.find_elements_by_class_name('highlight')
        # the above does not work for all cases
        
        # elem=driver.find_elements_by_xpath('//table[@class="jig-ncbigrid gene-tabular-rprt ui-ncbigrid"]/tbody/tr/td/div[2]/a')
        # the above form is deprecated; use form below
        try:
            elem=driver.find_element(by='xpath',value='//table[@class="jig-ncbigrid gene-tabular-rprt ui-ncbigrid"]/tbody/tr/td/div[2]/a')
        except:
            OutDict[name]['gene_id']='not found'
            OutDict[name]['ensembl_id']='not found'
            OutDict[name]['species']='not found'
            print(f'WARNING: {name} not found.')
            continue

        gene_id=elem.text

        # click the gene name
        elem.click()
        
        # locate search result; save search results as a list var
        elem=driver.find_element(by='id',value='summaryDl')
        s=elem.text.split('\n')
        
        ensembl_id=[item for item in s if 'nsembl' in item][0][-18:].split(' ')[0]
        species=s[s.index([item for item in s if 'rganism' in item][0])+1]

        #OutDict[name].append(gene_id, ensembl_id, organism)
        OutDict[name]['gene_id']=gene_id
        OutDict[name]['ensembl_id']=ensembl_id
        OutDict[name]['species']=species

        print(f'{name}, {OutDict[name]}')
        #driver.back()
        #driver.back()
    driver.quit()
    return OutDict


d2=convert(l,d, 'rhesus')
print(d2)

# output the string for pasting into R
#s=''
#for m in d:
#    if d[m]!=():
#        for marker in d[m]:
#            print(f'"{marker}",',end='')
