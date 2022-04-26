from bitarray import bitarray




class FingerPrint:

    def __init__(self, d_input):

        self.d_input = d_input


    def prepFP(self, l_bit_order = []):

        bits_out = bitarray()
        
        if l_bit_order != []:
            self.l_name_bit = l_bit_order
            
        else:
            self.l_name_bit = list(self.d_input.keys())

        # check if k are in it
        if "INPUT" in self.l_name_bit:
            self.l_name_bit.remove("INPUT")
        if "DTXSID" in self.l_name_bit:
            self.l_name_bit.remove("DTXSID")
        if "PREFERRED_NAME" in self.l_name_bit:
            self.l_name_bit.remove("PREFERRED_NAME")

        for name_bit in self.l_name_bit:
            bits_out.extend(self.d_input[name_bit])
        
        self.bits = bits_out

    
    def bitsToStrings(self, bitsIn=""):
        
        if bitsIn != "":
            stringbit = str(self.bits).split("'")[1]

        if not "bits" in self.__dict__:
            return ""
        else:
            stringbit = str(self.bits).split("'")[1]
        return stringbit


    
    def jaccardScore(self, bitsToCompare, exclude0=0):


        if len(self.bits) != len(bitsToCompare):
            print ("ERROR: no same length of bits - l.6 calculate")
            return "ERROR"

        #print("==check FP==")
        #print(self.bitsToStrings(self.bits))
        #print(self.bitsToStrings(bitsToCompare))

        i = 0
        identic = 0.0
        while i < len(self.bits):
            if self.bits[i] == bitsToCompare[i]:
                if exclude0 == 1:
                    if self.bits[i] != 0:
                        identic += 1.0
                else:    
                    identic += 1.0
            i += 1

        if identic != 0.0:
            score = identic / (len(self.bits) + len(bitsToCompare) - identic)
        else:
            score = 0.0

        return score


    def DICEScore(self, bitsToCompare, exclude0=0):

        if len(self.bits) != len(bitsToCompare):
            print ("ERROR: no same length of bits - l.6 calculate")
            return "ERROR"


        i = 0
        identic = 0.0
        while i < len(self.bits):
            if self.bits[i] == bitsToCompare[i]:
                if exclude0 == 1:
                    if self.bits[i] != 0:
                        identic += 1.0
                else:    
                    identic += 1.0
            i += 1

        if identic != 0.0:
            score = (2*identic) / (len(self.bits) + len(bitsToCompare))
        else:
            score = 0.0

        return score

