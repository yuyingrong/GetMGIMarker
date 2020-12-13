
# YR 20201130
# Updated: YR 20201213 minor bug fix
# Python 3.8.3
# this script scrapes gene ID info from: http://www.informatics.jax.org/batch/summary
# by modifying the list of input protein names,
# running the script will return and print their gene IDs on the console
# Note: this script runs very slowly

from selenium import webdriver
import time

d={}
# enter custom markers; below are my examples
l=['EMR1','Cd66b','Cd66b','Cd3','Cd4','Cd8','Cd45R','B220','Cd19','Cd22','Cd11c','Cd123','NKp46','Cd34','Cd11b','Cd66b','Cd41','Cd61','Cd9','Cd62P','Cd235a','Cd146','Cd106','Cd31','Cd62E','Cd326']


def convert(NamesList,MarkersDict):
    print('Note: This is a very slow program as it pauses 4 seconds in each run.\n')
    # opens a Chrome browser; must have Chrome web browser downloaded
    # specify path to the web browser
    browser=webdriver.Chrome('C:\Ming\Software\chromedriver_win32\chromedriver.exe')
    # opens the Jax webpage in the browser
    browser.get('http://www.informatics.jax.org/batch/summary')
    browser.implicitly_wait(10)
    
    for name in NamesList:
        # create an entry for the protein name in the dictionary
        MarkersDict[name]=()
        
        # locate the textarea, whose <id> is 'ids'
        elem=browser.find_element_by_id('ids')
        # clear the previous search term if any
        elem.clear()
        
        # type protein name into textarea
        elem.send_keys(name)
        # submit the form containing this <id>; i.e. search the protein ID
        elem.submit()

        # force Python to wait 4 sec so that the JS table is fully loaded
        # super inefficient and stupid step, but works for now
        time.sleep(4)
        
        # locate search result; save search results as a list var
        ColList=browser.find_elements_by_class_name('yui-dt-col-symbol')
        # var c contains >= 2 items, items.text = 'symbol' (the header), 'marker_name1', 'marker_name2', ...

        # check if the returned marker name is empty; "no associated gene" found by the website
        if ColList[-1].text=='':
            print(f'Searching {name} returned an empty marker name.  {name} will remain in the name-marker dictionary, but not appear in the string for pasting into R.')
        else:
            # append the marker name[s] into MarkersDict
            for marker in ColList[1:]:
                MarkersDict[name]+=(marker.text,)
        
        print(f'{name},{MarkersDict[name]}')
        browser.back()
    browser.quit()


convert(l,d)


# output the string for pasting into R
s=''
for m in d:
    if d[m]!=():
        for marker in d[m]:
            print(f'"{marker}",',end='')

# need to know if mult identical marker names are put in R
# e.g. both 'Cd45R' and 'B220' translate to 'Ptprc'
# will entering 'Ptprc' twice in R cause a problem for Seurat?

