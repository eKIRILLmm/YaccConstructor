module GLLAbstractParserBioTest

open System.IO
open QuickGraph
open NUnit.Framework
open AbstractAnalysis.Common
open Yard.Generators.GLL.AbstractParser
open Yard.Generators.Common.ASTGLL
open Yard.Generators.Common.ASTGLLFSA
open Yard.Generators.GLL.ParserCommon
open YC.API
open Yard.Frontends.YardFrontend
open Yard.Generators.GLL
open Yard.Core.Conversions.ExpandMeta

open System.Collections.Generic
open System.Linq

let dataDir = @"C:/projects/YC/YaccConstructor/tests/data/AbstractGLL_LabelledVert/"
let grammarsDir = @"C:/projects/YC/YaccConstructor/tests/GLL.AbstractParser.Simple.Tests/"

let getInputGraphVertLbl tokenizer inputFile startPos =    
    let edges = 
        File.ReadAllLines (dataDir + inputFile)
        |> Array.filter(fun x -> not (x = ""))
        |> Array.map (fun s -> let x = s.Split([|' '|])
                               x.[0], x.[1], x.[2])
    let edg (f : string) (t : string) (l : string) = 
        new TaggedEdge<_,_>(f, l, t)

    let g = new GraphLabelledVertex<_>(startPos, (fun x -> ((tokenizer x) |>int)))
    
    [|for (first, tag, last) in edges -> edg first tag last |]
    |> g.AddEdges
    |> ignore
    
    g 

let getParserSource grammarFile conv = 
    let fe = new YardFrontend()
    let gen = new GLL()
    generate (grammarsDir + grammarFile)
             fe gen 
             None
             conv
             [|""|]
             [] :?> ParserSourceGLL

let test grammarFile inputFile startPos nodesCount edgesCount termsCount ambiguityCount = 
    //printfn "%A" needChangeDirectory
    let conv = [new ExpandMeta()]
    let parser = getParserSource grammarFile conv
    let input  = getInputGraphVertLbl parser.StringToToken inputFile startPos
    let tree = buildAst parser input
    //printfn "%A" tree
    tree.AstToDot (dataDir + inputFile + ".dot")
    let n, e, t, amb = tree.CountCounters
    //printfn "%d %d %d %d" n e t amb
    Assert.AreEqual(nodesCount, n, sprintf "Nodes expected:%i, found:%i." nodesCount n)
    Assert.AreEqual(edgesCount, e, sprintf "Edges expected:%i, found:%i." edgesCount e)
    Assert.AreEqual(termsCount, t, sprintf "Terms expected:%i, found:%i." termsCount t) 
    Assert.AreEqual(ambiguityCount, amb, sprintf "Ambiguities expected:%i, found:%i." ambiguityCount amb)
    Assert.Pass()

[<TestFixture>]
type ``GLL abstract parser graph lbl vert tests``() =
    [<Test>]  
    member this._04_RightRecursionCheck () =
        test "RightRecursionCheck.yrd" 
             "RightRecursionCheck.txt"
             [|"P"|] 15 16 4 1

[<EntryPoint>]
let main argv = 
    System.Runtime.GCSettings.LatencyMode <- System.Runtime.GCLatencyMode.LowLatency
    let t = new ``GLL abstract parser graph lbl vert tests``() 
    test "RightRecursionCheck.yrd" 
         "RightRecursionCheck.txt"
         [|"P"|] 15 16 4 1
//    BioDataPreproc.preprocBioData()
//    BioDataPerformance.performTests()
    0