module BioDataPerformance

open VDS.RDF
open VDS.RDF.Parsing

open QuickGraph
open AbstractAnalysis.Common
open Yard.Generators.GLL.AbstractParser
open Yard.Generators.Common.ASTGLL
open Yard.Generators.GLL.ParserCommon

open System.IO
open System.Collections.Generic



let (|Prefix|_|) (p:string) (s:string) =
    if s.StartsWith(p) 
    then Some() 
    else None

let (|Equals|_|) x y = 
    if x = y 
    then Some() 
    else None

let getTokenFromTag tokenizer (tag:string) = 
    match tag with
    | Prefix "Protein_" () -> tokenizer "PROTEIN"
    | Prefix "Gene_" () -> tokenizer "GENE"
    | Prefix "Phenotype_" () -> tokenizer "PHENOTYPE"
    | Equals "interacts_with" () -> tokenizer "INTERACTS" 
    | Equals "belongs_to" () -> tokenizer "BELONGS" 
    | Equals "-belongs_to" () -> tokenizer "RBELONGS" 
    | Equals "codes_for" () -> tokenizer "CODESFOR"
    | Equals "-codes_for" () -> tokenizer "RCODESFOR"
    | Equals "refers_to" () -> tokenizer "REFERS"
    | Prefix "GO_" () -> tokenizer "GO"
    | Prefix "FamilyOrDomain_" () -> tokenizer "FAM_OR_DOM"
    | Equals "has_FamilyOrDomain" () -> tokenizer "HAS_FAM_OR_DOM"
    | _ -> tokenizer "OTHER"

let getEdgesVert file =
    let mutable count = 0
    let lines = File.ReadLines(file)
    [|
    for l in lines ->
//        if count = 0 || count = 1
//        then 
//            count <- count + 1
            let elems = l.Split('\t')
            elems.[0], elems.[1], elems.[2]
//        else
//            if count = 25
//            then count <- 0
//            else count <- count + 1
//            [||]
    |] 
//    |> Array.concat

let getEdges file startPos = 
    let vMap = new System.Collections.Generic.Dictionary<_,_>()
    let mutable startPosId = -1<positionInInput>
    let mutable idV = -1
    let getId v = 
        if vMap.ContainsKey v 
        then true, vMap.[v] 
        else (idV <- idV + 2; vMap.Add(v, idV); false, idV)
    let mutable count = 0
    let lines = File.ReadLines(file)
    [|
    for l in lines ->
//        if count = 0 || count = 1
//        then 
//            count <- count + 1
            let elems = l.Split('\t')
            let f = elems.[0]
            let lbl = elems.[1]
            let t = elems.[2]
            
            let getFId = getId f
            let getTId = getId t

            if f = startPos 
            then 
                match getFId with
                    | (_, id) -> startPosId <- id * 1<positionInInput>
            if t = startPos 
            then 
                match getTId with
                    | (_, id) -> startPosId <- id * 1<positionInInput>


            match getFId, getTId with
            | (false, fId), (false, tId) -> 
                [|fId, f, fId + 1;
                fId + 1, lbl, tId;
                tId, t, tId + 1;|]
                                            
            | (true, fId), (false, tId) ->
                [|fId + 1, lbl, tId;
                tId, t, tId + 1;|]

            | (false, fId), (true, tId) ->
                [|fId, f, fId + 1;
                fId + 1, lbl, tId;|]

            | (true, fId), (true, tId) -> 
                [|fId + 1, lbl, tId;|]
//        else
//            if count = 25
//            then count <- 0
//            else count <- count + 1
//            [||]

    |] |> Array.concat, startPosId


let getParseInputGraph file =
    
    let edges, startPos = getEdges file "Gene_2"
//    let allVs = edges |> Array.collect (fun (f,l,t) -> [|f * 1<positionInInput>; t * 1<positionInInput>|]) |> Set.ofArray |> Array.ofSeq
    let allVs = edges |> Array.choose (fun (f,l,t) -> 
        if getTokenFromTag id l = "GENE"
        then Some(f * 1<positionInInput>)
        else None) |> Set.ofArray |> Array.ofSeq
//    printfn "%A" startPos
    let graph = new SimpleInputGraph<_>(allVs, id)
    
    (*edges
    |> Array.collect (fun (f,l,t) -> [|new ParserEdge<_>(f, t, getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l)|])
    |> graph.AddVerticesAndEdgeRange
    |> ignore*)

    edges
    |> Array.collect (fun (f,l,t) -> 
//        if getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l <> (int) GLL.BioCFG.stringToToken.["OTHER"]
//        then [|new ParserEdge<_>(f, t, getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l)|]
//        else [||])
          [|new ParserEdge<_>(f, t, getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l)|])
    |> graph.AddVerticesAndEdgeRange
    |> ignore

    graph, graph.EdgeCount

let getParseInputGraphVert file =
    let edges = getEdgesVert file

    let allGenes = edges |> Array.choose (fun (f,l,t) -> 
        if l = "codes_for"
        then Some(f)
        else None)
    let startPos = allGenes |> Set.ofArray |> Array.ofSeq
//    let graph = new GraphLabelledVertex<int>([|getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) "Gene_348"|], [||], id)
    let graph = new GraphLabelledVertex<string>(startPos, [||], getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]))

    edges 
    |> Array.collect (fun (f,l,t) -> [|new TaggedEdge<_,_>( f, t, l)|])
    |> graph.AddEdges
    |> ignore

    graph, graph.EdgeCount
        
let processFile file =
    let g1, edges = 
        getParseInputGraph file 

    printfn "Number of edges: %A" edges
    let start = System.DateTime.Now
    

    let gss, sppf, _ = parse  GLL.BioCFG.parserSource g1 true
//    let roots = sppf.GetRoots gss 1<positionInInput>
    let nt = sppf.GetNonTermByName "s" GLL.BioCFG.parserSource
    let pathset = sppf.Iterate nt GLL.BioCFG.parserSource 1000000

    let a = sppf.GetTerminalNodes
    let b = sppf.TerminalNodes
    let gr = pathset.ToAdjacencyGraph
//    let root1 =
//        Yard.Generators.GLL.AbstractParser.getAllSPPFRoots GLL.BioCFG.parserSource g1
//    root1.[0].AstToDot "resultV.dot"
//    let root1 =
//        Yard.Generators.GLL.AbstractParser.getAllRangesForStartState GLL.GPPerf1.parserSource g1
//        |> Seq.length

    let time1 = (System.DateTime.Now - start).TotalMilliseconds
//    printfn "roots %A" root1.Length
    edges, time1//, root1

let performTests() =
    let allTriplesFile = @"..\..\..\data\BioData\result\allTriples.txt"
    let simpleInputFile = @"..\..\..\data\BioData\result\simpleInput.txt"
    let simpleInputFile1 = @"..\..\..\data\BioData\result\simpleInput1.txt"
    let simpleInputFile2 = @"..\..\..\data\BioData\result\simpleInput2.txt"
    processFile simpleInputFile1    
    |> printfn "%A"
    
    printfn "finished"
    System.Console.ReadKey() |> ignore