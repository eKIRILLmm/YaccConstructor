namespace AbstractLexer.Core

open System.Collections.Generic
open AbstractLexer.Common
open QuickGraph.Algorithms
open Microsoft.FSharp.Collections
open QuickGraph.Graphviz
open AbstractParsing.Common

[<Struct>]
type StateInfo<'a, 'b> =
    val StartV: int
    val AccumulatedString: ResizeArray<'a>
    val BackRefs: 'b
    new (startV, str, brs) = {StartV = startV; AccumulatedString = str; BackRefs = brs}

[<Struct>]
type State<'a, 'b> =
    val StateID: int
    val AcceptAction: int
    val Info: ResizeArray<StateInfo<'a, 'b>>
    new (stateId, acceptAction, info) = {StateID = stateId; AcceptAction = acceptAction; Info = info}

[<Sealed>]
type AsciiTables(trans: uint16[] array, accept: uint16[]) =
    let rec scanUntilSentinel(lexBuffer, state) =
        let sentinel = 255 * 256 + 255 
        // Return an endOfScan after consuming the input 
        let a = int accept.[state] 
        if a <> sentinel then () 
            //onAccept (lexBuffer,a)
        
        // read a character - end the scan if there are no further transitions 
        let inp = 1//int(lexBuffer.Buffer.[lexBuffer.BufferScanPos])
        let snew = int trans.[state].[inp] 
        if snew = sentinel then ()
            //lexBuffer.EndOfScan()
        else 
            //lexBuffer.BufferScanLength <- lexBuffer.BufferScanLength + 1;
            scanUntilSentinel(lexBuffer, snew)
            
    /// Interpret tables for an ascii lexer generated by fslex. 
    member tables.Interpret(initialState,lexBuffer) = 
        //startInterpret(lexBuffer)
        //scanUntilSentinel(lexBuffer, initialState)
        1

    static member Create(trans,accept) = new AsciiTables(trans,accept)

[<Sealed>]
type UnicodeTables(trans: uint16[] array, accept: uint16[]) = 
    let sentinel = 255 * 256 + 255 
    let numUnicodeCategories = 30 
    let numLowUnicodeChars = 128 
    let numSpecificUnicodeChars = (trans.[0].Length - 1 - numLowUnicodeChars - numUnicodeCategories)/2
    let lookupUnicodeCharacters (state,inp) = 
        let inpAsInt = int inp
        // Is it a fast ASCII character?
        if inpAsInt < numLowUnicodeChars then
            int trans.[state].[inpAsInt]
        else 
            // Search for a specific unicode character
            let baseForSpecificUnicodeChars = numLowUnicodeChars
            let rec loop i = 
                if i >= numSpecificUnicodeChars then 
                    // OK, if we failed then read the 'others' entry in the alphabet,
                    // which covers all Unicode characters not covered in other
                    // ways
                    let baseForUnicodeCategories = numLowUnicodeChars+numSpecificUnicodeChars*2
                    let unicodeCategory = System.Char.GetUnicodeCategory(inp)
                    //System.Console.WriteLine("inp = {0}, unicodeCategory = {1}", [| box inp; box unicodeCategory |]);
                    int trans.[state].[baseForUnicodeCategories + int32 unicodeCategory]
                else 
                    // This is the specific unicode character
                    let c = char (int trans.[state].[baseForSpecificUnicodeChars+i*2])
                    //System.Console.WriteLine("c = {0}, inp = {1}, i = {2}", [| box c; box inp; box i |]);
                    // OK, have we found the entry for a specific unicode character?
                    if c = inp
                    then int trans.[state].[baseForSpecificUnicodeChars+i*2+1]
                    else loop(i+1)
                
            loop 0    
        
    let scanUntilSentinel inp (state:State<_,_>) =
        // Return an endOfScan after consuming the input 
        let a = int accept.[state.StateID]
        let onAccept = if a <> sentinel then a else state.AcceptAction
        // Find the new state
        let snew = lookupUnicodeCharacters (state.StateID,inp)
        snew = sentinel, onAccept, snew

    let tokenize actions (states:array<_>) edgesSeq lastVId printG =
        let edgesSeq = edgesSeq |> Array.ofSeq
        let add (edg:AEdge<_,_>) (newStt:State<_,_>) =
            match states.[edg.Target]
                  |> ResizeArray.tryFind(fun (x:State<_,_>) -> x.AcceptAction = newStt.AcceptAction && x.StateID = newStt.StateID)
                with
            | Some x ->
                newStt.Info
                |> ResizeArray.iter(
                    fun i -> 
                        if x.Info.Exists(fun j -> j.StartV = i.StartV
                                                  && i.AccumulatedString.Count = j.AccumulatedString.Count
                                                  && ResizeArray.forall2 (=) i.AccumulatedString j.AccumulatedString) 
                            |> not
                        then x.Info.Add i)
            | None -> states.[edg.Target].Add newStt
 
        let mkNewString (edg:AEdge<_,_>) (stt:State<_,_>) =
            let ch = edg.Label.Value
            if stt.Info.Count > 0
            then
                stt.Info
                |> ResizeArray.map 
                    (
                        fun i -> 
                            new StateInfo<_,_>(i.StartV
                            , ResizeArray.concat [i.AccumulatedString; ResizeArray.singleton ch]
                            , ResizeArray.concat [i.BackRefs; ResizeArray.singleton edg.BackRef])
                    )
            else 
                new StateInfo<_,_>(edg.Source, ResizeArray.singleton ch, ResizeArray.singleton edg.BackRef)
                |> ResizeArray.singleton

        let processEdg (edg:AEdge<_,_>) stt reduced =
            let acc = new ResizeArray<_>(10) 
            match edg.Label with
            | Some x ->
                let rec go stt =
                    let reduce, onAccept, news = scanUntilSentinel x stt
                    if reduce
                    then
                        for i in stt.Info do 
                            actions onAccept (new string(i.AccumulatedString |> Array.ofSeq)) (i.BackRefs |> Array.ofSeq)
                            |> fun x -> 
                                if not !reduced then acc.Add(new ParserEdge<_>(i.StartV,edg.Source, x))
                                reduced := true
                        let newStt = new State<_,_>(0,-1,new ResizeArray<_>())                            
                        go newStt
                    else 
                        let acc = mkNewString edg stt
                        let newStt = new State<_,_>(news,onAccept,acc)
                        add edg newStt
                go stt      
            | None -> add edg stt

            acc 

        let res_edg_seq = 
            seq{
                for (edgs:array<AEdge<_,_>>) in edgesSeq do
                    for stt in states.[edgs.[0].Source] do
                        let reduced = ref false
                        for edg in edgs do
                            yield! processEdg edg stt reduced
                }
        seq{
            yield! res_edg_seq
            for x in states.[lastVId] do
                for i in x.Info do                        
                    let x = actions (int accept.[x.StateID]) (new string(i.AccumulatedString.ToArray())) (i.BackRefs.ToArray()) 
                    yield (new ParserEdge<_>(i.StartV,lastVId, x))
        }
        

    let inputGraph actions (inG:LexerInputGraph<_>) printG = 
        let g = new LexerInnerGraph<_>(inG)
        let sorted = g.TopologicalSort() |> Array.ofSeq
        let states = Array.init ((Array.max sorted)+1) (fun _ -> new ResizeArray<_>())
        let startState = new State<_,_>(0,-1, ResizeArray.singleton (new StateInfo<_,_>(0,new ResizeArray<_>(), new ResizeArray<_>())))
        states.[g.StartVertex] <- ResizeArray.singleton startState
        let edgesSeq = seq{ for v in sorted do
                              yield g.OutEdges v |> Array.ofSeq
                                  
                           }
                       |> Seq.filter (fun x -> x.Length > 0)
        let newEdgs = tokenize actions states edgesSeq sorted.[sorted.Length-1] printG |> Array.ofSeq
        let res = new ParserInputGraph<_>()
        let r = newEdgs
        res.AddVerticesAndEdgeRange r
        |> ignore
        res
                          
    // Each row for the Unicode table has format 
    //      128 entries for ASCII characters
    //      A variable number of 2*UInt16 entries for SpecificUnicodeChars 
    //      30 entries, one for each UnicodeCategory
    //      1 entry for EOF


    member tables.Tokenize(actions,g, ?printG) = inputGraph actions g (match printG with Some f -> f | _ -> fun x y -> ())
    static member Create(trans,accept) = new UnicodeTables(trans,accept)

//    let query = "select \" f , " + (if x < 2 then "x" else "y") + "from z"
//    let DB = new MS_DB("")
//    let res = DB.exec query
//    res |> Seq.iter (printfn "%A")