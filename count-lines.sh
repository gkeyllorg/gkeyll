# count lines in libary and app
wc -l core/apps/*.h core/apps/*.c core/creg/*.c core/unit/*.c core/zero/*.h core/zero/*.c \
   moments/apps/*.h moments/apps/*.c moments/creg/*.c moments/luareg/*.lua moments/unit/*.c moments/zero/*.h moments/zero/*.c \
   vlasov/apps/*.h vlasov/apps/*.c vlasov/creg/*.c vlasov/luareg/*.lua vlasov/unit/*.c vlasov/zero/*.h vlasov/zero/*.c \
   gyrokinetic/apps/*.h gyrokinetic/apps/*.c gyrokinetic/creg/*.c gyrokinetic/luareg/*.lua gyrokinetic/unit/*.c gyrokinetic/zero/*.h gyrokinetic/zero/*.c \
   pkpm/apps/*.h pkpm/apps/*.c pkpm/creg/*.c pkpm/luareg/*.lua pkpm/unit/*.c pkpm/zero/*.h pkpm/zero/*.c \
    | sort -n
