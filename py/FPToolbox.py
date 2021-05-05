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
            print(self.d_input[name_bit])
            bits_out.extend(self.d_input[name_bit])
        
        self.bits = bits_out

    
    def bitsToStrings(self):
        
        if not "bits" in self.__dict__:
            return ""

        stringbit = str(self.bits).split("'")[1]
        print(stringbit)
        return stringbit


    
    def jaccardIndex(self, bitsToCompare):


        if len(self.bits) != len(bitsToCompare):
            print ("ERROR: no same length of bits - l.6 calculate")
            return "ERROR"

        i = 0
        identic = 0.0
        while i < len(self.bits):
            if self.bits[i] == bitsToCompare[i]:
                identic += 1.0
            i += 1

        if identic != 0.0:
            score = identic / (len(self.bits) + len(bitsToCompare) - identic)
        else:
            score = 0.0

        return score

