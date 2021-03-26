/*
A collection of helper methods to make an HTML summary.
TODO: Implement similar functions for latex
*/
TString accEffiTableConverter(std::vector<std::vector<TString>> accEffiVector)
{
    TString htmlTable = "<table class='customTable' style='width:100%'>";
    int numberOfRows = (int)accEffiVector.size();
    for (int iRow = 0; iRow < numberOfRows; iRow++)
    {
        htmlTable.Append("<tr>\n");
        std::vector<TString> rowVector = accEffiVector[iRow];
        int numberOfCols = (int)rowVector.size();
        for (int iCol = 0; iCol < numberOfCols; iCol++)
        {
            if (iRow == 0)
            {
                htmlTable.Append(Form("<th>%s</th>\n", rowVector[iCol].Data()));
            }
            else
            {
                htmlTable.Append(Form("<td>%s</td>\n", rowVector[iCol].Data()));
            }
        }
        htmlTable.Append("</tr>\n");
    }
    htmlTable.Append("</table>");
    return htmlTable;
}

TString takeTemplateAndReplace(TString templatePath, std::map<TString, TString> dictionary)
{
    TString original = gSystem->GetFromPipe(Form("cat %s", templatePath.Data()));
    TString output = original;
    std::map<TString, TString>::iterator it;

    for (it = dictionary.begin(); it != dictionary.end(); it++)
    {
        output = output.ReplaceAll(it->first, it->second);
    }
    return output;
}
void MakeSummaryPerCentrality(int centMin, int centMax, int numberOfPerformedIterations, Double_t largestDiffAccEffiToPrevious, std::vector<std::vector<TString>> accEffiVectorVsRap, std::vector<std::vector<TString>> accEffiVectorVsPt)
{
    ofstream outputFile(Form("./Cent-%dto%d/Summary.html", centMin, centMax), std::ofstream::out);
    outputFile << Form("<h2 id='cent_%d-%d'>Centrality %d-%d</h2>", centMin, centMax, centMin, centMax) << endl;
    outputFile << takeTemplateAndReplace("HTML/Snipets/SummarySection.html", {{"#NUMEBROFITERATIONS", Form("%d", numberOfPerformedIterations)}, {"#SECONDLASTITERATION", Form("%d", numberOfPerformedIterations - 1)}, {"#LARGESTDIFFTOPREVIOUS", Form("%2.2f", largestDiffAccEffiToPrevious)}}) << endl;
    outputFile << takeTemplateAndReplace("HTML/Snipets/SummaryFigures.html", {{"#CENTDIR", Form("Cent-%dto%d", centMin, centMax)}}) << endl;
    outputFile << "<h5> Acceptance&#215;efficiency values for various iterations as a function of pt: </ h5> " << endl;
    outputFile << accEffiTableConverter(accEffiVectorVsPt) << endl;
    outputFile << "<h5> Acceptance&#215;efficiency values for various iterations as a function of rapidity: </ h5> " << endl;
    outputFile << accEffiTableConverter(accEffiVectorVsRap) << endl;
    outputFile << "<h3> Plots per iteration step : </ h3> " << endl;
    outputFile << takeTemplateAndReplace("HTML/Snipets/Step-0_Plots.html", {{"#CENTDIR", Form("Cent-%dto%d", centMin, centMax)}}) << endl;

    for (int iterationStep = 1; iterationStep <= numberOfPerformedIterations; iterationStep++)
    {
        outputFile << takeTemplateAndReplace("HTML/Snipets/Step-i_Plots.html", {{"#CENTDIR", Form("Cent-%dto%d", centMin, centMax)}, {"#ITERATIONSTEP", Form("%d", iterationStep)}, {"#PREVIOUSITERATIONSTEP", Form("%d", iterationStep - 1)}}) << endl;
    }
    outputFile.close();
}

void CollectCentSummaries(vector<vector<int>> centralityBins)
{
    TDatime dateNow;
    ofstream outputFile(Form("./HTML/Summary.html"), std::ofstream::out);
    int numberOfCentralityBins = (int)centralityBins.size();
    TString mainDiv = Form("<h1>Analysis on %s</h1>\n", dateNow.AsString());
    TString sideNav = "<p>Centralities</p>";

    for (int iCent = 0; iCent < numberOfCentralityBins; iCent++)
    {
        mainDiv.Append(gSystem->GetFromPipe(Form("cat ./Cent-%dto%d/Summary.html", centralityBins[iCent][0], centralityBins[iCent][1])));
        sideNav.Append(Form("<a href='#cent_%d-%d'>%d-%d%% </a>\n", centralityBins[iCent][0], centralityBins[iCent][1], centralityBins[iCent][0], centralityBins[iCent][1]));
    }

    TString summary = takeTemplateAndReplace("./HTML/Snipets/Document.html", {{"#ALLCENTSUMMARIES", mainDiv}, {"#SIDENAV", sideNav}});
    outputFile << summary << endl;
    outputFile.close();

    gSystem->Exec(Form("open ./HTML/Summary.html"));
}