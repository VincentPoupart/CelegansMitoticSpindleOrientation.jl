function SpindleOrientationCellShape(path)
    ### Make a foilder where results will be saved
    PathResults = mkpath(joinpath(path, "results"))

    #### fetching the raw data from the path for mitotic and interphase cells
    filelist = ListFiles(path)
    wrlfiles = sort(WrlfileList(filelist))
    DTC = FetchDTC(filelist)
    DPaxis = FetchDPaxis(filelist)
    CellShape_df = FetchImarisStats(filelist, DTC, DPaxis)

    ### Add angles betwween mean long axis and long axis for each time point
    CellShape_df = AnglesMeanCtoC(CellShape_df)

    ### Create Mitosis DataFrame
    Mitosis_dict = CreateCelloutput(filelist)

    ### Scoring Mitotic Stages (needs user inputs)
    scoring_df = ClickStepsMitosis(Mitosis_dict)

    ### Create Interphase_df by keeping only stable and constant cells (see filtration parameters)
    Interphase_df = FilteringCellShape(CellShape_df)

    ### Add surface normal data to Interphase
    Interphase_df = NormalVector2InterphaseDF!(Interphase_df, wrlfiles)

    ### Add surface normal data to Mitosis
    Mitosis_dict = NormalVector2MitosisDict!(Mitosis_dict, wrlfiles)
    
    #### Dict to DataFrame
    Mitosis_df = CelloutputDict2DF(Mitosis_dict)

    ### Add cell shape track and mitosis steps to Mitosis_df
    Mitosis_df = FindCellShapeTracksAndMitoticSteps(Mitosis_df, CellShape_df, scoring_df, DPaxis)

    ### Remove missing Data
    Mitosis_df = dropmissing(Mitosis_df, [:AnaphaseOnset, :MeanCellVolume])

    ### Filtering Mitosis_df
    Mitosis_df = FilteringCelloutput(Mitosis_df)

    ### Remove Mitotic Data to keep only interphase cells
    Interphase_df = RemoveMitoticTrackID(Mitosis_df, Interphase_df)

    ### Add surface normal data to Interphase_df
    Interphase_df = NormalVector2InterphaseDF!(Interphase_df, wrlfiles)

    ### Add DP data
    Mitosis_df = AddDPaxisMitosis!(Mitosis_df, DPaxis)
    Interphase_df = AddDPaxisInterphase!(Interphase_df, DPaxis)
    Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset] = convert.(Float64, Mitosis_df[!, :MeanAngleSpindleLongAxisAnaphaseOnset])
    Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset] = convert.(Float64, Mitosis_df[!, :MeanAngleSpindleRachisAnaphaseOnset])
    Mitosis_df = dropmissing(Mitosis_df, :MeanAngleSpindleDPaxis)

    ### Add projected axis data to Interphase_df
    Interphase_df = AddProjectAxesInterphaseDF!(Interphase_df)

    ### Add projected axis data to Mitosis_df
    Mitosis_df = AddProjectAxesMitosisDF!(Mitosis_df)

    ### Clean Mitosis_df by removing Missing values
    Mitosis_df = dropmissing(Mitosis_df, [:MeanAngleSpindleRachisAnaphaseOnset, :MeanAngleSpindleLongAxisAnaphaseOnset, :MeanOrientationCtoDP, :MeanAngleLongAxisRachis, :MeanAngleProjectLongestAxisDP, :MeanAngleProjectSpindleProjectLongAxis, :MeanAngleProjectSpindleProjectShortAxis])

   
    ### Save Results
    CSV.write(joinpath(PathResults, "InterphaseDF.csv"), Interphase_df)
    CSV.write(joinpath(PathResults, "MitosisDF.csv"), Mitosis_df)

end