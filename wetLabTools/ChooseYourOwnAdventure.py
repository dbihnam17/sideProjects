##This program is meant to automate the choose your own adventure sheet
#Cell counts are based off ideal amounts required for sequencing
#This program is meant to be used POST-SORT, so PRE-SORT scRNA isn't accounted for

#CONSTANTS

IDEAL_SCRNA = 1.0
IDEAL_PDNA = 0.5
IDEAL_RNA = 2.0
IDEAL_ATACLO = 0.5
IDEAL_ATACHI = 1.0
IDEAL_PSCDNAHI = 1.0
IDEAL_PSCDNALO = 0.5
IDEAL_NDNA = 1.0
IDEAL_NSCDNA = 1.0
P_TOTAL = 5.0

#Sets

projectSet = {'PROJECT','PROJECT1','PROJECT 1','PROJECT2','PROJECT 2'}
req24Set = {'REQ24','REQ 24'}
req25Set = {'REQ25','REQ 25'}
batchSet = {'PROJECT 1, PROJECT 2, Req24, Req25'}

def main() :
    
    #Collect batch name
        
    batch = str(input(('Enter sample batch name: ')))
    dash40()
    batch = batch.upper() #Eliminate case sensitivity in batch names

    #Calculate sequencing aliquots with functions

    if batch in projectSet or batch in req25Set:
        posCell, negCell = cellCount()
        pDNA, rna, atac, scDNA, nDNA, check, warning = projectFunction(posCell, negCell)
        outMessageFunction(batch, pDNA, rna, atac, scDNA, nDNA, 0, check, warning)
    elif batch in req24Set:
        posCell, negCell = cellCount()
        pDNA, rna, atac, nscRNA, check, warning = req24Function(posCell, negCell)
        outMessageFunction(batch, pDNA, rna, atac, 0, 0, nscRNA, check, warning)
    else: #Input validation for batch name
        dash63()
        print('Failed to provide correct or supported batch name')
        print('Supported batches: %s' %(batchSet))
        dash63()

    #Function to collect cell counts

def cellCount():
    print('Please enter cell count as so:\n5,500,000 cells should be input as 5.5') #Example input
    dash40()
    posCell = float(input('Enter double-positive cell count: '))
    negCell = float(input('Enter negative cell count: '))
    return posCell, negCell

#project function

def projectFunction(posCell, negCell):
    check = False
    warning = ''
    while True:
        if 0 <= posCell < 3:
            check = True
            pDNA, rna, atac, warning, scDNA = checkLow(posCell, check)
        elif 3.0 <= posCell < 3.5:
            pDNA = IDEAL_PDNA
            rna = (posCell - IDEAL_PDNA - IDEAL_ATACLO)
            atac = IDEAL_ATACLO
            scDNA = 0.0
        elif 3.5 <= posCell < 4.0:
            pDNA = IDEAL_PDNA
            rna = (posCell - IDEAL_PDNA - IDEAL_ATACHI)
            atac = IDEAL_ATACHI
            scDNA = 0.0
        elif 4.0 <= posCell < 5.0:
            pDNA = IDEAL_PDNA
            rna = (posCell - IDEAL_PDNA - IDEAL_ATACHI - IDEAL_PSCDNALO)
            atac = IDEAL_ATACHI
            scDNA = IDEAL_PSCDNALO
        elif 5.0 <= posCell:
            pDNA = (IDEAL_PDNA + ((posCell - P_TOTAL) / 3))
            rna = (IDEAL_RNA + (2 * (posCell - P_TOTAL) / 3))
            atac = IDEAL_ATACHI
            scDNA = IDEAL_PSCDNAHI
        else: #Input validation for positive cells
            print('%.2f is not a valid input for double-positive cell count.' %(posCell))
            posCell = float(input('Enter double-positive cell count: '))
            continue
        break
    while True: #Input validation for negative cells
        if 0 <= negCell:
            nDNA = negCell
        else:
            print('%.2f is not a valid input for negative cell count.' %(negCell))
            negCell = float(input('Enter negative cell count: '))
            continue
        break
    return pDNA, rna, atac, scDNA, nDNA, check, warning

#Req24 function

def req24Function(posCell, negCell):
    check = False
    warning = ''
    while True:
        if 0 <= posCell < 3:
            check = True
            pDNA, rna, atac, warning, scDNA = checkLow(posCell, check)
        elif 3.0 <= posCell < 3.5:
            pDNA = IDEAL_PDNA
            rna = (posCell - IDEAL_PDNA - IDEAL_ATACLO)
            atac = IDEAL_ATACLO
        elif 3.5 <= posCell < 5.0:
            pDNA = IDEAL_PDNA
            rna = (posCell - IDEAL_PDNA - IDEAL_ATACHI)
            atac = IDEAL_ATACHI
        elif 5.0 <= posCell:
            pDNA = (IDEAL_PDNA + ((posCell - P_TOTAL) / 3))
            rna = (IDEAL_RNA + (2 * (posCell - P_TOTAL) / 3))
            atac = IDEAL_ATACHI
        else: #Input validation for positive cells
            print('%.2f is not a valid input for double-positive cell count.' %(posCell))
            posCell = float(input('Enter double-positive cell count: '))
            continue
        break
    while True: #Input validation for negative cells
        if 0 <= negCell:
            nscRNA = negCell
        else:
            print('%.2f is not a valid input for negative cell count.' %(negCell))
            negCell = float(input('Enter negative cell count: '))
            continue
        break
    return pDNA, rna, atac, nscRNA, check, warning

#Output function

def outMessageFunction(batch, pDNA, rna, atac, scDNA, nDNA, nscRNA, check, warning):
    dash40()
    print()
    print('Bulk DNA: %.2fM cells' %(pDNA))
    if warning == '*Order extra vial for RNA' :
        print('Bulk RNA: %.2fM cells*' %(rna))
    else :
        print('Bulk RNA: %.2fM cells' %(rna))
    print('ATAC: %.2fM cells' %(atac))
    if batch in projectSet or batch in req25Set:
        print('scDNA: %.2fM cells' %(scDNA))
        print('Negative DNA: %.2fM cells' %(nDNA))
        print()
    elif batch in req24Set:
        print('Negative scRNA: %.2fM cells' %(nscRNA))
        print()
    if check:
        print(warning)
        print()
    dash40()

#Function to handle low ++ cell count cases and check
#If we have more vials to pull from
    
def checkLow(posCell, check) :
    exVialSet = {'Y', 'N', 'YES', 'NO'}
    if check :
        dash40()
        exVial = input('Can you order another vial? (Y/N): ').upper()
        inExVialSet = False
        while not inExVialSet :
            if exVial in exVialSet :
                inExVialSet = True
                if exVial == 'Y' or exVial == 'YES' :
                    warning = '*Order extra vial for RNA'
                    scDNA = 0
                    if 0 <= posCell < 0.8 :
                        pDNA = posCell
                        rna = 0
                        atac = 0
                    elif 0.8 <= posCell < 1 :
                        pDNA = (posCell - IDEAL_ATACLO)
                        rna = 0
                        atac = IDEAL_ATACLO
                    elif 1 <= posCell < 3 :
                        pDNA = IDEAL_PDNA
                        rna = (posCell - IDEAL_PDNA - IDEAL_ATACLO)
                        atac = IDEAL_ATACLO
                    else :
                        return print('ERROR: Invalid cell count in checkLow')
                    return pDNA, rna, atac, warning, scDNA
                else :
                    warning = 'Consult with Esteban'
                    scDNA = 0
                    if 0 <= posCell < 0.8 :
                        pDNA = posCell
                        rna = 0
                        atac = 0
                    elif 0.8 <= posCell < 1 :
                        pDNA = (posCell - IDEAL_ATACLO)
                        rna = 0
                        atac = IDEAL_ATACLO
                    elif 1 <= posCell < 3 :
                        pDNA = IDEAL_PDNA
                        rna = (posCell - IDEAL_PDNA - IDEAL_ATACLO)
                        atac = IDEAL_ATACLO
                    else :
                        return print('ERROR: Invalid cell count in checkLow')
                    return pDNA, rna, atac, warning, scDNA
            else:
                print('Invalid input.')
                exVial = input('Can you order another vial? (Y/N): ').upper()
    else :
        return print('ERROR: No check entry in checkLow')

#Function to print 40 dashes

def dash40() :
    return print('-' * 40)

#Function to print 63 dashes

def dash63() :
    return print('-' * 63)
    
main()