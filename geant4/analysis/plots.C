#include <vector>
#include <iostream>

void plots(const char* filename)
{
	TFile* tf = new TFile(filename);
	TTree* Tracks = tf->Get<TTree>("events");
	TH2D* hist_bar_layer = new TH2D("hist_bar_layer", "hits per bar per layer", 25, 0., 24., 17, 0., 16.);

	std::vector<int>* layerID{};
	Tracks->SetBranchAddress("layerID", &layerID);
	std::vector<int>* barID{};
	Tracks->SetBranchAddress("barID", &barID);
	std::vector<int>* EcalEdep{};
	Tracks->SetBranchAddress("EcalEdep", &EcalEdep);

	int number_of_entries = Tracks->GetEntries();
	for (int n = 0; n < number_of_entries; n++) {
		Tracks->GetEntry(n);
		int size = layerID->size();
		for (int i = 0; i < size-1; i++) {
			if ((*EcalEdep)[i] > 0.) {
				hist_bar_layer->Fill((*barID)[i], (*layerID)[i]);
			}
		}
	}
	TCanvas* c_bar_layer = new TCanvas("c_bar_layer", "c_bar_layer", 1200, 600);
	hist_bar_layer->SetXTitle("bar #");
	hist_bar_layer->SetYTitle("layer #");
	hist_bar_layer->SetZTitle("hits");
	hist_bar_layer->GetZaxis()->SetTitle("counts");
	// hist_bar_layer->SetContour(1000);
	hist_bar_layer->Draw("colz");
	c_bar_layer->SaveAs("bar_layer.pdf");
	c_bar_layer->SaveAs("bar_layer.C");
	tf->Close();
}
