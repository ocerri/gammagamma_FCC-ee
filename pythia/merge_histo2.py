import ROOT as rt
import glob

h_dsdeta_f = rt.TH1F()
h_dsdpt_f = rt.TH1F()
h_had_per_ev_f = rt.TH1F()
h_W_f = rt.TH1F()
h_id_part_f = rt.TH1F()

file_list = glob.glob('_root/MSTP14-10/*.root')



n_file=0
    file = rt.TFile(file_name)
    ## if n_file==0:
    ##     h_dsdeta_f = file.Get('hdsigmadeta')
    ##     h_dsdpt_f = file.Get('hdsigmadpT')
    ##     h_had_per_ev_f = file.Get('h_hadrons_per_event')
    ##     h_W_f = file.Get('hW')
    ##     h_id_part_f = file.Get('h_id_particle')
    ##     print(type(h_dsdeta_f))
    ## else:
    print(type(h_dsdeta_f))
    h_dsdeta_f.(file.Get('hdsigmadeta'))
    h_dsdpt_f.Merge(file.Get('hdsigmadpT'))
    h_had_per_ev_f.Merge(file.Get('h_hadrons_per_event'))
    h_W_f.Merge(file.Get('hW'))
    h_id_part_f.Merge(file.Get('h_id_particle'))
    n_file+=1
    file.Close()
    print(type(h_dsdeta_f))

h_dsdeta_f.Scale(1/n_file)
h_dsdpt_f.Scale(1/n_file)

out_file = rt.TFile('_root/MSTP14-10_final.root')
out_file.cd()
h_dsdeta_f.Write()
h_dsdpt_f.Write()
h_had_per_ev_f.Write()
h_W_f.Write()
h_id_part_f.Write()
out_file.Close()
