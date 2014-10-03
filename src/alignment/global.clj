(ns alignment.global
  "Implementation of a global sequence alignment algorithm for strings.

   This implementation is mostly educational, and can be optimized significantly:

   - the clojure/java itself can be optimized further.

   - there are algorithmic improvements in the literature

   This particular implementation stems from work originally undertaken by
   Needleman and Wunsch [*] and as taught by the textbook: HANDBOOK OF
   COMPUTATIONAL MOLECULAR BIOLOGY [**]

   [*, **] Citations in resource directory.")


   ;; - TODO: Improve documentation
   ;; - TODO: Optimize code
   ;; - TODO: Optimize algorithm logic.
   ;; - TODO: Add citations


;; Example Comparision Strings
(def ^:private a (into [] "agcttg"))
(def ^:private b (into [] "aggctga"))

;; Scoring Rules
(def ^:private conserved 2)
(def ^:private mismatch -1)
(def ^:private gap -1)


(defn- table
  "Build a 'table' (vector of vectors) of x and y dimension"
  [x y]
  (vec
   (for [i (range y)]
     (into [] (repeat x nil)))))

(defn- scond-y0
  "Transform a 'table' to satisfy the global alignment table
   starting conditions for table values where the y coordinate = 0."
  [align-table]
  (assoc align-table 0 (into []
                             (map #(* mismatch %)
                                  (-> align-table first count range)))))


(defn- scond-x0
  "Transform a 'table' to satisfy the global alignment table
   starting conditions for table values where the x coordinate
   = 0."
  ([align-table]
     (if (empty? align-table)
       align-table
       (scond-x0 align-table (dec (count align-table)))))
  ([align-table decrement]
     (if (= decrement 0)
       align-table
       (scond-x0
        (assoc-in align-table [decrement 0] (* mismatch decrement))
        (dec decrement)))))


(defn- arbitrary-ga-table
  "Create a 'global alignment table' configured with the starting
   conditions for the algorithm. Dimensions are arbitrary values.

   S[0,0] = 0
   S[0,y] = y * mismatch_penalty
   S[x,0] = x * mismatch_penalty"
  [x y]
  (-> (table x y) scond-x0 scond-y0))


(defn init-ga-table
  "Create an initialized 'gloabal alignment table' that satisfies
   the starting conditions for the global alignment algorithm.
   Accepts the sequences, and sizes table accordingly. The return
   value should be passed into the algorithm's fn chain"
  [seqa seqb]
  (let [genseq (fn [s] (into [] s))
        a (genseq seqa)
        b (genseq seqb)]
    (-> (table (inc (count a)) (inc (count b)))
        scond-x0
        scond-y0)))
