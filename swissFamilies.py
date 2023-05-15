## from https://www.swisssidechain.ch/data/family_table.pdf
Phenylalanine_derivatives=['TFP3', 'PTR', 'F2F', 'NAL', 'OTYR', 'TFP4', 'PTR2', 'DBY', 'NAO1', 'MTY', '4BF', 'PBF', 'TYQ', 'NAO2', 'MPH2', 'PHI', 'BIF', 'TYI', 'ANTH', 'APD', 'CNP2', 'DAH', 'PF5', 'PYR2', '4PH', '3CF', '3MY', 'BB8', 'PYR3', 'TBP4', '4CF', 'IYR', 'OMX', 'PYR4', 'CPH2', 'HOX', 'TY2', 'OMY', 'QU32', 'FCL', '0A1', 'YOF', 'HPE', 'QU33', '200', 'APM', 'NIY', 'DIPH', 'QU34', 'FPH2', '0BN', 'MP34', 'BCS', 'QU35', 'FPH3', '4HMP', 'CP24', 'STYA', 'QU36', 'PFF', 'DMP3', 'CP34', 'KYN', 'HQA', 'TFP2', 'PPN', 'WFP', 'ALN', 'QX32']

Phenylglycine_derivatives=['004', 'D4P', '', 'CPG2', 'CPG3', 'CPG4', 'FPG2', 'FPG3', '', 'FPG4', 'TFG2', 'TFG3', 'TFG4', '3FG', 'CHP']

Tryptophan_derivatives=['CTE', '0AF', '6CW', 'TRX', 'BTR', 'HRP', '32T', '4HT', 'BTH3', 'TRO', 'HTR', 'TTQ', '4IN', 'MTR6', 'MTR5', 'MOT5', 'FT6', 'FTR', '4FW']

Methionine_derivatives=['ME0', 'ESC', '4CY', '2FM', 'SME', 'MHO', 'OMT']

ALA_derivatives=['2AG', 'ABA', 'ACZ', 'ADAM', 'AHP', 'ALC', 'BIU', 'BUG', 'C2N', 'CHG', 'CPA3', 'DILE', 'FLA', 'FVAL', 'HLEU', 'I2M', 'IGL', 'IIL', 'LEF', 'LVG', 'NLE', 'NVA', 'OBF', 'TFLE']
ASN_GLN_derivatives=['AHB', 'MEN', 'HGA', 'LMQ', 'YCM', 'GHG', 'MEQ']
ARG_derivatives=['GDPR', 'GBUT', 'HRG', 'ARO', 'AGM', 'GGB', 'CIR', 'THIC']
HIS_derivatives=['2HF', '2HF1', '2HF2', 'PYZ1', 'TRZ4', 'TEZA', 'TZA4', 'THA3', 'TIH', 'FUA2']
other_derivatives=['AS2', 'AZDA', 'CAN', 'CSA', 'LED', 'LVN', 'M2S', 'OAS', 'OLT', 'ONL', 'THG2', 'THG3']
LYS_derivatives=['DPP', 'DAB', 'ORN', 'HHK', 'SLZ']
SER_THR_derivatives=['AA4', 'ALO', 'CTH', 'DDZ', 'HIL4', 'HL2', 'HLU', 'HSER', 'HVA', 'LDO', 'ILX', 'OCY', 'SEP', 'SEP2', 'TH6', 'TPO', 'TPO2', 'VAH']
CYS_derivatives=['HCS', 'LE1']
ASP_GLU_derivatives=['BHD', '2AS', 'DMK', 'FGL', '3GL', 'LME', 'MEG', 'SYM', 'FGA4', 'UN1', '2NP', '26P', 'CCS', '6CL', 'GME']

coarse = {'F':Phenylalanine_derivatives + Phenylglycine_derivatives,
          'W':Tryptophan_derivatives,
          'M':Methionine_derivatives,
          'A':ALA_derivatives,
          'R':ARG_derivatives,
          'H':HIS_derivatives,
          'L':LYS_derivatives,
          'S':SER_THR_derivatives,
          'C':CYS_derivatives,
          'Q':ASN_GLN_derivatives,
          'E':ASP_GLU_derivatives,
          'o':other_derivatives,
          }
          
f = open('data/rotamers/swiss.lib', 'r')
lines = f.readlines()
f.close()

nbRot = 0
nbType = 0
notfound = []
for line in lines:
    if line.startswith('rotamer'):
        w = line.split()
        coarse_type=None
        nbRot += 1
        for k,v in coarse.items():
            if w[1] in v:
                coarse_type = k
                print('%4s %s'% (w[1], k))
                nbType += 1
        if coarse_type is None:
            notfound.append(w[1])

assert len(notfound)==0
assert nbType==nbRot

f = open('newswisslib.lib', 'w')

for line in lines:
    if line.startswith('rotamer'):
        w = line.split()
        for k,v in coarse.items():
            if w[1] in v:
                coarse_type = k
                break
        f.write('%s %s %s %s %s\n'%(w[0], w[1], coarse_type, w[2], w[3]))
    else:
        f.write('%s'%line)
                
f.close()
