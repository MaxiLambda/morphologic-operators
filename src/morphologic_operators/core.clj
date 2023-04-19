(ns morphologic-operators.core
    (:import (clojure.lang PersistentVector)))

(defn uniform [val height length]
    "length and height have to be positive integers"
    (let [row (vec (take length (repeat val)))]
        (vec (take height (repeat row)))))

(defn ones [height length]
    (uniform 1 height length))

(defn zeros [height length]
    (uniform 0 height length))

(defn pretty [matrix]
    (doseq [row matrix]
        (println row)))

(defn vectorize [col]
    (vec (map vec col)))

(defn to-matrix [length col]
    (vectorize (partition length col)))

(defn get-at
    ([matrix y x]
     (-> matrix
         (get y)
         (get x)))
    ([matrix yx]
     (get-at matrix (first yx) (second yx))))

(defn is-valid? [max-y max-x yx]
    (let [[y x] yx]
        (and
            (< y max-y)
            (> y -1)
            (< x max-x)
            (> x -1))))

(defn errosion-check [image yx]
    "checks if the image shall be eroded based on the given indices"
    (not (every? #(= 1 %) (map #(get-at image %) yx))))

(defmulti invert (fn [arg] (class arg)))
(defmethod invert Number [val]
    (bit-xor 1 val))

(defmethod invert PersistentVector [val]
    (vec (map invert val)))
(defn errosion
    "1 -> 0
    image and mask are both matrices containing only 0 and 1
    mask has to be:
    size mask <= size image
    odd length mask
    odd height mask"
    ([image mask]
     (let [image-max-x (count (nth image 0))
           image-max-y (count image)
           mask-max-x (count (nth mask 0))
           mask-max-y (count mask)
           mask-center-x (/ (- mask-max-x 1) 2)
           mask-center-y (/ (- mask-max-y 1) 2)]
         (to-matrix image-max-x
                    (for [y (range image-max-y)
                          x (range image-max-x)
                          :let [y'x' (for [y' (range mask-max-y)
                                           x' (range mask-max-x)]
                                         [y' x'])
                                ;coordinates in the mask with 1
                                y'x'-ones (filter #(= 1 (get-at mask %)) y'x')
                                ;transform to offsets in image
                                y'x'-offset (map #(map - % [mask-center-y mask-center-x]) y'x'-ones)
                                y'x'-image (map #(map + % [y x]) y'x'-offset)
                                ;look only at points in the image with valid coordinates
                                y'x'-valid (filter #(is-valid? image-max-y image-max-x %) y'x'-image)]]
                        (if (or (= 0 (get-at image y x)) (errosion-check image y'x'-valid)) 0 1)))))
    ([image] (errosion image (ones 3 3))))

(defn diletation-check [image yx] "checks if the image has a 1 at any given coordinate"
    (some? (some #(= 1 %) (map #(get-at image %) yx))))

(defn diletation
    "0 -> 1
    image and mask are both matrices containing only 0 and 1
    mask has to be:
    size mask <= size image
    odd length mask
    odd height mask"
    ([image mask]
     (let [image-max-x (count (nth image 0))
           image-max-y (count image)
           mask-max-x (count (nth mask 0))
           mask-max-y (count mask)
           mask-center-x (/ (- mask-max-x 1) 2)
           mask-center-y (/ (- mask-max-y 1) 2)]
         (to-matrix image-max-x
                    (for [y (range image-max-y)
                          x (range image-max-x)
                          :let [y'x' (for [y' (range mask-max-y)
                                           x' (range mask-max-x)]
                                         [y' x'])
                                ;coordinates in the mask with 1
                                y'x'-ones (filter #(= 1 (get-at mask %)) y'x')
                                ;transform to offsets in image
                                y'x'-offset (map #(map - % [mask-center-y mask-center-x]) y'x'-ones)
                                y'x'-image (map #(map + % [y x]) y'x'-offset)
                                ;look only at points in the image with valid coordinates
                                y'x'-valid (filter #(is-valid? image-max-y image-max-x %) y'x'-image)]]
                        (if (or (= 1 (get-at image y x)) (diletation-check image y'x'-valid)) 1 0)))))
    ([image] (diletation image (ones 3 3))))

(defn opening "removes clusters and noise"
    ([image mask]
     (-> image
         (errosion mask)
         (diletation mask)))
    ([image] (opening image (ones 3 3))))

(defn closing "removes holes"
    ([image mask]
     (-> image
         (diletation mask)
         (errosion mask)))
    ([image] (closing image (ones 3 3))))