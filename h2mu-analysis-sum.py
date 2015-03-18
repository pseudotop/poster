import ROOT
ROOT.gSystem.Load("libFWCoreFWLite.so");
ROOT.gSystem.Load("libDataFormatsFWLite.so");
ROOT.AutoLibraryLoader.enable()
from DataFormats.FWLite import Events, Handle
import os,sys

def looseJet(jet):
    out = False;
    if (jet.neutralHadronEnergy() + jet.HFHadronEnergy() ) / jet.energy() < 0.99 \
        and jet.neutralEmEnergyFraction() < 0.99 \
        and jet.numberOfDaughters() > 1 \
        and jet.chargedHadronEnergyFraction() > 0 or abs(jet.eta()) > 2.4 \
        and jet.chargedMultiplicity() > 0 or abs(jet.eta()) > 2.4 \
        and jet.chargedEmEnergyFraction() < 0.99 or abs(jet.eta()) > 2.4: 
        out = True
    return out

cut = [0,0,0,0,0,0]
def isTightMuon(muon, vtx):
    cut[0] += 1
    if not muon.isPFMuon():
        cut[1] += 1
        if not muon.isGlobalMuon():
            cut[2] += 1
            return False
    muID = (muon.numberOfMatchedStations() > 1)
    if muID: cut[3] += 1
    if ( muon.innerTrack().isNull() or muon.muonBestTrack().isNull() ) : return False
    hits = muon.innerTrack().hitPattern().trackerLayersWithMeasurement() > 5 and muon.innerTrack().hitPattern().numberOfValidPixelHits() > 0
    #ip = abs(muon.muonBestTrack().dxy(vtx.position())) < 0.2 and abs(muon.muonBestTrack().dz(vtx.position())) < 0.5
    #if ip: cut[4] += 1
    tight = muID #and ip
    if tight: cut[5] +=1
    return tight
Phase_target = ["Phase1_PU_140_1000fb_1", "Phase1_PU_50", "Phase2_PU_140"]

print("=" * 50)
print("1. %s" % Phase_target[0])
print("2. %s" % Phase_target[1])
print("3. %s" % Phase_target[2])
print("=" * 50)

num = input("which program do you want to use? : ")


#Events("DIR/h2mu_ggh----.root"), DIR list : Phase1_PU_140_1000fb_1  Phase1_PU_50  Phase2_PU_140  
#events = Events("/pnfs/user/h2mu_TP/Phase2_PU_140/h2mu_ggh_M125GeV_14TeV_Phase2_PU_140_RECO_10.root")
path = "/pnfs/user/h2mu_TP/%s/" % Phase_target[num-1]
f_list = os.listdir(path)
target = []
 
for x in f_list:
    if 'root' in x:
        tmp2 = x.split("_")[1]
        if tmp2 == "ggh":
            target.append(x)
print target

out_root_tot = ROOT.TFile(Phase_target[num-1]+".root", "RECREATE")
out_num_event_tot = open(Phase_target[num-1]+"_num_of_events.txt","w")

if not os.path.isdir("/home/lsupertopl/h2mu_ggh"):
    os.mkdir("/home/lsupertopl/h2mu_ggh")
os.chdir("/home/lsupertopl/h2mu_ggh")

out_tr = ROOT.TTree("TotTree","Sum of the Higgs Mass")

# histogram
hist_tot = ROOT.TH1F("hist", "Higgs Mass Distribution", 200, 0., 200.)
hist_tot.GetXaxis().SetTitle("GeV")
hist_tot.GetYaxis().SetTitle("# of Higgs")

hist2 = ROOT.TH1F("hist2", "# of muons per each event Distribution", 4, 0, 4)
mu_num = ["0", "1", "2", "3", "4"]
hist2.GetXaxis().SetTitle("# of muons per an event")
hist2.GetYaxis().SetTitle("# of events happening")
hist2.SetFillColor(38)
sbl = 0
while sbl < 5:
    sbl += 1
    hist2.GetXaxis().SetBinLabel(sbl, mu_num[sbl-1])

cut_muon_pt = 30.
cut_muon_eta = 2.1
cut_muon_iso = 10.12
cut_jet_pt = 30.
cut_jet_eta = 4.7

sum_tot = [0,0,0]
sum_iev = 0
for x in target:
    path_tot = path + x
    print("=" * 50) 
    print path_tot
    print("=" * 50) 
    events = Events(path_tot)

    div_dir = path_tot.split('/')
    out_f = div_dir[-1] 
    out_root = ROOT.TFile(out_f,"RECREATE")
    out_num_event = open(out_f[:-5]+"_num_of_events.txt","w")

    # histogram
    hist = ROOT.TH1F("hist", "Higgs Mass Distribution", 200, 0., 200.)
    hist.GetXaxis().SetTitle("GeV")
    hist.GetYaxis().SetTitle("# of Higgs")


    jetsLabel, jets = "ak5PFJets", Handle("std::vector<reco::PFJet>")
    muonsLabel, muons = "muons", Handle("std::vector<reco::Muon>")
    vertexsLabel, vertexs = "offlinePrimaryVertices", Handle("std::vector<reco::Vertex>")
    nmuons = 0
    njets = 0
    npass = 0
    for iev,event in enumerate(events):
        print iev
        selectedmuons=[]
        selectedjets=[]

        event.getByLabel(vertexsLabel, vertexs)
        vert = vertexs.product()[0]
        if not vert:
            print "no vertax"
        event.getByLabel(muonsLabel,muons)
        #print len(muons.product())
        a = [0,0,0,0,0]
        for g,m in enumerate(muons.product()):
            #if m.isPFMuon() and ( m.isGlobalMuon() or m.isTrackerMuon()):
            a[0] += 1
            if not isTightMuon(m, vert):
                continue
            a[1] += 1
            if m.pt() < cut_muon_pt:
                continue
            a[2] += 1
            if abs(m.eta()) > cut_muon_eta:
                continue
            a[3] += 1
            pfIsolationR04 = m.pfIsolationR04().sumChargedHadronPt+m.pfIsolationR04().sumPhotonEt+m.pfIsolationR04().sumNeutralHadronEt
            if pfIsolationR04/m.pt() > cut_muon_iso:
                continue
            a[4] += 1
            selectedmuons.append(m)
            nmuons +=1
        print a
        # the number of muons per event#  
        ev_mu = a[4]
        hist2.Fill(mu_num[ev_mu],1)
        ################################
        print cut
        event.getByLabel(jetsLabel,jets)
        for g,j in enumerate(jets.product()):
            if not looseJet(j):
                continue
            if j.pt() < cut_jet_pt:
                continue
            if abs(j.eta()) > cut_jet_eta:
                continue
            selectedjets.append(j)
            njets += 1
        if len(selectedmuons) < 2:
            continue

        for m in selectedmuons:
            print "muon pt", m.pt(), "muon eta", m.eta(), "muon charge", m.charge()
            
        if selectedmuons[0].charge()*selectedmuons[1].charge() == 1:
            continue
        mu1 = selectedmuons[0]
        mu2 = selectedmuons[1]
    
        higgs = ROOT.TLorentzVector(mu1.px(), mu1.py(), mu1.pz(), mu1.energy()) + ROOT.TLorentzVector(mu2.px(), mu2.py(), mu2.pz(), mu2.energy())

        print higgs.M()
        h_M = higgs.M()
    
        hist.Fill(h_M)
        hist_tot.Fill(h_M)
        npass +=1
        
    print "nmuons",nmuons
    print "njets",njets
    print "npass",npass

    sum_iev += iev+1

    sum_tot[0] += nmuons
    sum_tot[1] += njets
    sum_tot[2] += npass

    sum_cut = []
    sum_cut += cut

    hist.Draw()
    out_root.Write()
    out_root.Close()    

print sum_cut
print sum_tot
hist_tot.Draw()
#hist2.LabelsDeflate()
#hist2.LabelsOption("u","X")
hist2.Draw()
out_root_tot.Write()
out_root_tot.Close()

out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n # of events : %d\n total # of muons : %d\n total # of jets : %d\n total # of pass : %d\n" % (sum_iev,sum_tot[0],sum_tot[1],sum_tot[2]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n /* muID = (muon.numberOfMatchedStations() > 1), ip = abs(muon.muonBestTrack().dxy(vtx.position())) < 0.2 and abs(muon.muonBestTrack().dz(vtx.position())) < 0.5, tight = muID and ip */\n ")
out_num_event_tot.write(" in isTightMuon func. total tight-muon processes\n # of nocut : %d\n # of non-ParticleFlow : %d\n # of non-GlobalMuon : %d\n" % (cut[0],cut[1],cut[2]))
out_num_event_tot.write(" # of muID : %d\n # ip : %d\n # of tight : %d\n" % (cut[3],cut[4],cut[5]))
out_num_event_tot.write("=" * 50)
out_num_event_tot.write("\n")
out_num_event_tot.close()
