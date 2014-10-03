(ns alignment.global
  "Implementation of a global sequence alignment algorithm for strings.

   This implementation is mostly educational, and can be optimized significantly:

   - the clojure/java itself can be optimized further.

   - there are algorithmic improvements in the literature

   This particular implementation stems from work originally undertaken by
   Needleman and Wunsch [*] and as taught by the textbook: HANDBOOK OF
   COMPUTATIONAL MOLECULAR BIOLOGY [**]

   This is an original Clojure implementation. All of the code is an interpretation
   of the general concepts described in the textbook. This code is the original
   work of the author (bitfl1pper). Obviously, the code is merely a derivative of the
   logic described in the textbook.

   The maths of the algorithm can be found on pages # #2 #3.

   I plan on licensing this code to be as free as possible. But I don't really
   know much about licensing, so that's on my todo list.

   For now, assume Eclispse Public License. This code must retain citations.
   Most importantly: Needleman / Wunsch et al., the Handbook [...], and this
   author.

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


(defn paircomp
  "I am unwilling to document this function at the moment.

   Create the pair/nopair truth table."
  ([seqa seqb]
     (let [comp-char (first seqa)
           comp (map #(= comp-char %) seqb)
           accumulate (into [] comp)]
       (paircomp (rest seqa) seqb accumulate)))
  ([seqa seqb accumulate]
     (if (empty? seqa)
       (vec (map #(vec %) (partition (count seqb) accumulate)))
       (let [comp-char (first seqa)
             comp (map #(= comp-char %) seqb)
             carry (into accumulate comp)]
         (paircomp (rest seqa) seqb carry)))))

(defn score-lookup [x y align-table]
  (get-in align-table [x y]))

(defn nw-score-lookup [x y align-table]
  (score-lookup (dec x) (dec y) align-table))

(defn n-score-lookup [x y align-table]
  (score-lookup x (dec y) align-table))

(defn w-score-lookup [x y align-table]
  (score-lookup (dec x) y align-table))

(defn match?score [x y align-table pair-table]
  (cond
   (= x 0)
   (get-in align-table [y x])
   (= y 0)
   (get-in align-table [y x])
   (true? (get-in pair-table [(dec y) (dec x)]))
   conserved
   :else
   mismatch))


(defn diag?score [x y align-table pair-table]
  (+ (nw-score-lookup x y align-table)
     (match?score x y align-table pair-table)))

(defn n?score [x y align-table]
  (+ (n-score-lookup x y align-table)
     mismatch))

(defn w?score [x y align-table]
  (+ (w-score-lookup x y align-table)
     mismatch))

(defn score [x y align-table pair-table]
  (con (= x 0) (score-lookup x y align-table)
       (= y 0) (score-lookup x y align-table)
       :else
       (max (diag?score x y align-table pair-table)
            (n?score x y align-table)
            (w?score x y align-table))))

(defn table-dimsension [table]
  {:x (count (first table))
   :y (count table)})

(defn has-score? [x y align-table]
  (not (nil? (score-lookup x y align-table))))

(defn update-atable [align-table pair-table x y]
  (assoc-in align-table [x y] (score x y align-table pair-table)))
