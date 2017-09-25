module BioDataPerformance

open VDS.RDF
open VDS.RDF.Parsing

open QuickGraph
open AbstractAnalysis.Common
open Yard.Generators.GLL.AbstractParser
open Yard.Generators.Common.ASTGLL
open Yard.Generators.GLL.ParserCommon
open YC.GLL.SPPF

open System.IO
open System.Collections.Generic

let PrintToDotVert (graph: AdjacencyGraph<_,TaggedEdge<_,_>>) name (tagToString: _ -> string) (*(tokenToString : 'token -> string) (numToToken : int -> 'token)*) = 
    use out = new System.IO.StreamWriter (name : string)
    out.WriteLine("digraph AST {")
    out.WriteLine "rankdir=LR"
//        for i=0 to graph.VertexCount-1 do
//            out.Write (i.ToString() + "; ")
//        out.WriteLine()
    for i in graph.Vertices do
        let edges = graph.OutEdges i
        for e in edges do
            let tokenName = e.Tag |> tagToString 
            out.WriteLine (e.Source.ToString() + " -> " + e.Target.ToString() + "[label=\"" + tokenName + "\"]")
    out.WriteLine("}")
    out.Close() 

let SubgraphToGraphVert (subgraph: AdjacencyGraph<_,ParserEdge<_>>) (input: IParserInput) =
    let graphVert = new AdjacencyGraph<_,_>()
    for i in subgraph.Vertices do
        if i % 2 = 0
        then
            let vFrom = input.PositionToString i
            let edges1 = subgraph.OutEdges i
            for e1 in edges1 do
                let edges2 = subgraph.OutEdges e1.Target
                for e2 in edges2 do
                    let vTo = input.PositionToString e2.Target 
                    graphVert.AddVerticesAndEdge (new TaggedEdge<_,_>(vFrom, vTo, e2.Tag)) |> ignore
    graphVert

let PrintToDot (graph: AdjacencyGraph<_,ParserEdge<_>>) name (tagToString: _ -> string) (*(tokenToString : 'token -> string) (numToToken : int -> 'token)*) = 
    use out = new System.IO.StreamWriter (name : string)
    out.WriteLine("digraph AST {")
    out.WriteLine "rankdir=LR"
//        for i=0 to graph.VertexCount-1 do
//            out.Write (i.ToString() + "; ")
//        out.WriteLine()
    for i in graph.Vertices do
        let edges = graph.OutEdges i
        for e in edges do
            let tokenName = e.Tag |> tagToString 
            out.WriteLine (e.Source.ToString() + " -> " + e.Target.ToString() + "[label=\"" + tokenName + "\"]")
    out.WriteLine("}")
    out.Close() 

let private fst (f, _, _) = f

let private snd (_, s, _) = s

let private trd (_, _, t) = t

let private SPPFToSubgraph (sppf : SPPF) (ps : ParserSourceGLL) =
    let tagToLabel x = ps.IntToString.Item (x |> int)
    let edges = GetTerminals sppf |> Seq.map(fun x -> new ParserEdge<_>(snd x, trd x, (fst x |> tagToLabel)))
    let subgraph = new AdjacencyGraph<int, ParserEdge<_>>()
    subgraph.AddVerticesAndEdgeRange(edges) |> ignore
    subgraph

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
    let lines = File.ReadLines(file)
    [|
    for l in lines ->
        let elems = l.Split('\t')
        elems.[0], elems.[1], elems.[2]
    |]

let getEdges file = 
    let vMap = new System.Collections.Generic.Dictionary<_,_>()
    let mutable idV = -2
    let getId v = 
        if vMap.ContainsKey v 
        then true, vMap.[v] 
        else (idV <- idV + 2; vMap.Add(v, idV); false, idV)
    let mutable count = 0
    let lines = File.ReadLines(file)
    [|
    for l in lines ->
        let elems = l.Split('\t')
        let f = elems.[0]
        let lbl = elems.[1]
        let t = elems.[2]
            
        let getFId = getId f
        let getTId = getId t

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
    |] |> Array.concat


let getParseInputGraph file =
    
    let edges = getEdges file

    let allVs = edges |> Array.collect (fun (f,l,t) -> [|f * 1<positionInInput>; t * 1<positionInInput>|]) |> Set.ofArray |> Array.ofSeq
    let allGenes = edges |> Array.choose (fun (f,l,t) -> 
        if getTokenFromTag id l = "GENE"
        then Some(f * 1<positionInInput>)
        else None)
    let genes = allGenes |> Set.ofArray |> Array.ofSeq 

    let graph = new SimpleInputGraph<_>(allGenes, id)
    
    (*edges
    |> Array.collect (fun (f,l,t) -> [|new ParserEdge<_>(f, t, getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l)|])
    |> graph.AddVerticesAndEdgeRange
    |> ignore*)

    edges
    |> Array.collect (fun (f,l,t) -> 
        if getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l <> (int) GLL.BioCFG.stringToToken.["OTHER"]
        then [|new ParserEdge<_>(f, t, getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) l)|]
        else [||])
    |> graph.AddVerticesAndEdgeRange
    |> ignore

    graph, graph.EdgeCount

let getParseInputGraphVert file =
    let edges = getEdgesVert file
    let allGenes = edges |> Array.choose (fun (f,l,t) -> 
        if getTokenFromTag id f = "GENE"
        then Some(f)
        else None)
    let genes = allGenes |> Set.ofArray |> Array.ofSeq 
    let halfGenes = genes 
                    |> Array.indexed 
                    |> Array.choose (fun (i, x) -> 
                            if i % 99 = 0
                            then Some(x)
                            else None)
        
    let graph = new GraphLabelledVertex<string>(halfGenes, halfGenes, (fun t -> getTokenFromTag (fun x -> (int) GLL.BioCFG.stringToToken.[x]) t))

    edges 
    |> Array.collect (fun (f,l,t) -> [|new TaggedEdge<_,_>(f,t,l)|])
    |> graph.AddEdges
    |> ignore

    graph, graph.EdgeCount
        
let processFile file =
    let g1, edges = 
        getParseInputGraphVert file 

    printfn "%A" edges

    let start = System.DateTime.Now
   
    let _,sppf,_ = parse GLL.BioCFG.parserSource g1 true
    
    let time1 = (System.DateTime.Now - start).TotalMilliseconds

    PrintToDotVert g1 "inputGraph.dot" id

    let subgraph = SPPFToSubgraph sppf GLL.BioCFG.parserSource

    let subgraphVert = SubgraphToGraphVert subgraph g1 
    PrintToDotVert subgraphVert "subgraph.dot" id

    edges, time1

let performTests() =
    let allTriplesFile = @"..\..\..\data\BioData\result\allTriples.txt"
    let simpleInputFile = @"..\..\..\data\BioData\result\simpleInput.txt"
    let genes20 = @"..\..\..\data\BioData\result\20genes.txt"
    let genes50 = @"..\..\..\data\BioData\result\50genes.txt"
    let genes100 = @"..\..\..\data\BioData\result\100genes.txt"
    let genes300 = @"..\..\..\data\BioData\result\300genes.txt"
    let genes500 = @"..\..\..\data\BioData\result\500genes.txt"
    let genes700 = @"..\..\..\data\BioData\result\700genes.txt"
    let allGenes = @"..\..\..\data\BioData\result\allGenes.txt"
    let genes1000 = @"..\..\..\data\BioData\result\1000genes.txt"
    let genes2000 = @"..\..\..\data\BioData\result\2000genes.txt"
    let genes3000 = @"..\..\..\data\BioData\result\3000genes.txt"
    let genes5000 = @"..\..\..\data\BioData\result\5000genes.txt"
    let genes10000 = @"..\..\..\data\BioData\result\10000genes.txt"
    processFile genes300
    |> printfn "%A"
    
    printfn "finished"
    System.Console.ReadKey() |> ignore